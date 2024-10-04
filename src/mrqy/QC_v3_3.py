#############################################################################################################################
#############################################################################################################################
####################                                                                                      ###################
####################               MRQy - Quality control and evaluation tool for MRI data                ###################
####################                                                                                      ###################
#############################################################################################################################
#############################################################################################################################

'''
MRQy (v3) is a quality control and evaluation tool for MRI data.
This program first uses AI-assited brain extraction techniques to compute segmentation masks of the brain. You have the option
of using either FSL-BET (fast and a bit accurate) or HD-BET (slower but much more accurate).
Then, this mask is used to separate the foreground image (brain) from the background image (skull and others).
Then, several quality control metrics are extracted from the image metadata and other are calculated thanks to the foreground
and the background image.
Finally, AI data dimensional reduction algorithms (UMAP and TSNE) are used to help find outliers in the data and to visualize
the data in a lower dimension.
All the metrics and info are saved in a .csv (.tsv) file. To visualize the data, open the index.html file in a navigator 
(UserInterface/index.html). And then open the .tsv file. You will find tables, charts and graphs to  navigate through the 
dataset.
To function, the input folder must contain NIFTI images (.nii.gz). It doesn't function with DICOM images (.dcm) yet.
'''

#############################################################################################################################
##############################                        Import libraries                         ##############################
#############################################################################################################################

import argparse                                             # Handles command line arguments
import datetime                                             # Manipulate time values
import math                                                 # Mathematical functions
import matplotlib.cm as cm                                  # Interactive colormaps
import matplotlib.pyplot as plt                             # Interactive plots
from medpy.io import load                                   # For NIFTI image processing
import nibabel as nib                                       # Read and write access to common neuroimaging file formats
import numpy as np                                          # Mathematical operations over arrays
import os                                                   # Interaction with the operating system
import pandas as pd                                         # Data analysis
from scipy import ndimage                                   # Multidimensional image processing
from scipy.cluster.vq import whiten                         # Clustering algorithms
from scipy.ndimage import rotate                            # Rotation of arrays
from scipy.signal import convolve2d as conv2                # Signal processing, convolution
from skimage import exposure as ex
from skimage.filters import median                          # Image local median
from skimage.morphology import square                       # Square shaped footprint of an image
from sklearn.manifold import TSNE                           # t-SNE clustering, dimensionality reduction
import subprocess                                           # use command-line in python
import time                                                 # Manipulate time values
import umap                                                 # UMAP clustering, dimensionality reduction
import warnings                                             # Warning subsystem
warnings.filterwarnings("ignore")                           # Remove all warnings like conversion thumbnails

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

# Initialisation
nfiledone = 0
csv_report = None
first = True
headers = []


def print_msg_box(msg, indent=1, width=None, title=None):
    
    # Split the message into lines
    lines = msg.split('\n')
    
    # Define the space with indentation
    space = " " * indent
    
    # Determine the width based on the maximum line length if width is not specified
    if not width:
        width = max(map(len, lines))
    
    # Create the top border of the box
    box = f'╔{"═" * (width + indent * 2)}╗\n'  
    
    # Add title section if title exists
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n' 
    
    # Add message lines to the box
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    
    # Add the bottom border of the box
    box += f'╚{"═" * (width + indent * 2)}╝' 
    
    # Print the box
    print(box)


#############################################################################################################################
##############################                           Definitions                           ##############################
#############################################################################################################################


#################### File pre-processing ####################
#############################################################


def file_name(root):
    
    # Gathering relevant file paths based on extensions. MRQy supports .dcm and .nii.gz files
    files = [os.path.join(dirpath, filename) for dirpath, _, filenames in os.walk(root)
             for filename in filenames
             if filename.endswith('.nii.gz')]
    
    # Separating files based on their extensions
    nifti_images = [i for i in files if i.endswith('.nii.gz')]

    # Extracting subject identifiers from the different files
    nifti_names = [os.path.basename(scan)[:os.path.basename(scan).index('.')] for scan in nifti_images]

    # Displaying the total number of identified subjects / or images
    print('The number of images is {}'.format(len(nifti_names)))
    
    # Returning various lists containing file paths
    return files, nifti_names, nifti_images


#################### NIFTI files processing ####################
################################################################


def brain_extraction(scan, output_masks):
    
    # Make the path for the output masks
    output_mask = os.path.join(output_masks, os.path.basename(scan))
    
    # FSL-BET
    subprocess.run(['bet', scan, output_mask, '-f', '0.5', '-g', '0', '-m', '-n'], check=True)
    
    '''
    # HD-BET   
    subprocess.run(['hd-bet', '-i', scan, '-o', output_mask, '-tta', '0', '-pp', '0', '-b', '0', '-mode', 'fast'], check=True)
    '''
    mask_path = output_mask.replace(".nii.gz", "_mask.nii.gz")

    return mask_path


def volume_nifti_masks(mask_path):
    
    # Transform mask into 3D volume
    mask_data, mask_header = load(mask_path)
    mask_shape = np.shape(mask_data)
    if mask_shape[0] < mask_shape[1] and mask_shape[0] < mask_shape[2]:         # Sagittal
        masks = [mask_data[i, :, :] for i in range(mask_shape[0])]
    elif mask_shape[1] < mask_shape[0] and mask_shape[1] < mask_shape[2]:       # Coronal
        masks = [mask_data[:, i, :] for i in range(mask_shape[1])]               
    elif mask_shape[2] < mask_shape[0] and mask_shape[2] < mask_shape[1]:       # Axial
        masks = [mask_data[:, :, i] for i in range(mask_shape[2])]
    else:                                                                           # 3D
        masks = [mask_data[:, :, i] for i in range(mask_shape[2])]
        
    mask_name = os.path.basename(mask_path).split('.')[0]
    
    return masks, mask_name, mask_header

'''
def saveThumbnails_nifti_masks(volume_masks, output_directory):
    
    os.makedirs(output_directory + os.sep + volume_masks[1])
    
    for i in range(len(volume_masks[0])):
        plt.imsave(output_directory + os.sep + volume_masks[1] + os.sep + volume_masks[1] + '(%d).png' % int(i+1), rotate(volume_masks[0][i], 90), cmap = cm.Greys_r)
    
    print('The number of %d masks are saved to %s' % (len(volume_masks[0]), output_directory + os.sep + volume_masks[1]))
'''

def volume_nifti(scan, name):
    
    # Loading image data and header
    image_data, image_header = load(scan)
    
    # Get image dimensions
    image_shape = np.shape(image_data)
    
    # Extract the 2D images from the 3D image data. The correct 3D image volume is created with the extracted 2D nifti slices
    if image_shape[0] < image_shape[1] and image_shape[0] < image_shape[2]:         # Sagittal
        images = [image_data[i, :, :] for i in range(image_shape[0])]
    elif image_shape[1] < image_shape[0] and image_shape[1] < image_shape[2]:       # Coronal
        images = [image_data[:, i, :] for i in range(image_shape[1])]               
    elif image_shape[2] < image_shape[0] and image_shape[2] < image_shape[1]:       # Axial
        images = [image_data[:, :, i] for i in range(image_shape[2])]
    else:                                                                           # 3D
        images = [image_data[:, :, i] for i in range(image_shape[2])]
    
    # Return the images, their name and image_header
    return images, name, image_header


def saveThumbnails_nifti(volume, output_directory):
    
    # Create a directory for images     # In "UserInterface/Data/output_folder_name/image_name_folder"
    os.makedirs(output_directory + os.sep + volume[1])
        
    # Save images as thumbnails, with a rotation
    for i in range(len(volume[0])):
        # Save the images. Because of the precedent image processing operations, the image was rotated 90° in anti-clockwise. So we apply a 90° clockwise rotation
        # In "UserInterface/Data/output_folder_name/image_name_folder/image_name(i).png"
        plt.imsave(output_directory + os.sep + volume[1] + os.sep + volume[1] + '(%d).png' % int(i+1), rotate(volume[0][i],90), cmap = cm.Greys_r)
    
    # print
    print('The number of %d images are saved to %s' % (len(volume[0]), output_directory + os.sep + volume[1]))    


class BaseVolume_nifti(dict):

    def __init__(self, output_directory, volume, ol, scan, sample_size, ch_flag):
        
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        
        # Add information to the output list
        self.addToPrintList("Image", volume[1], volume, ol, 170)
        self["outdir"] = output_directory
        self.addToPrintList("Name of Images", os.listdir(output_directory + os.sep + volume[1]), volume, ol, 100)   # Here are the images shown in "UserInterface/index.html"
        self.addToPrintList("VRX", format(volume[2].get_voxel_spacing()[0], '.2f'), volume, ol, 1)
        self.addToPrintList("VRY", format(volume[2].get_voxel_spacing()[1], '.2f'), volume, ol, 2)
        self.addToPrintList("VRZ", format(volume[2].get_voxel_spacing()[2], '.2f'), volume, ol, 3)
        self.addToPrintList("ROWS", np.shape(volume[0])[1], volume, ol, 4)
        self.addToPrintList("COLS", np.shape(volume[0])[2], volume, ol, 5)
        self["os_handle"] = volume[0]
        self.addToPrintList("NUM", len(volume[0]), volume, ol, 6)
        self.addToPrintList("ORIENTATION", orientation(scan), volume, ol, 7) 
        self.addToPrintList("MEAN", vol(volume, volume_masks, sample_size, "Mean", output_directory, ch_flag), volume, ol, 8)
        self.addToPrintList("RNG", vol(volume, volume_masks, sample_size, "Range", output_directory, ch_flag), volume, ol, 9)
        self.addToPrintList("VAR", vol(volume, volume_masks, sample_size, "Variance", output_directory, ch_flag), volume, ol, 10)
        self.addToPrintList("CV", vol(volume, volume_masks, sample_size, "CV", output_directory, ch_flag), volume, ol, 11)
        self.addToPrintList("CPP", vol(volume, volume_masks, sample_size, "CPP", output_directory, ch_flag), volume, ol, 12)
        self.addToPrintList("PSNR", vol(volume, volume_masks, sample_size, "PSNR", output_directory, ch_flag), volume, ol, 13)
        self.addToPrintList("SNR1", vol(volume, volume_masks, sample_size, "SNR1", output_directory, ch_flag), volume, ol, 14)
        self.addToPrintList("SNR2", vol(volume, volume_masks, sample_size, "SNR2", output_directory, ch_flag), volume, ol, 15)
        self.addToPrintList("SNR3", vol(volume, volume_masks, sample_size, "SNR3", output_directory, ch_flag), volume, ol, 16)
        self.addToPrintList("SNR4", vol(volume, volume_masks, sample_size, "SNR4", output_directory, ch_flag), volume, ol, 17)
        self.addToPrintList("CNR", vol(volume, volume_masks, sample_size, "CNR", output_directory, ch_flag), volume, ol, 18)
        self.addToPrintList("CVP", vol(volume, volume_masks, sample_size, "CVP", output_directory, ch_flag), volume, ol, 19)
        self.addToPrintList("CJV", vol(volume, volume_masks, sample_size, "CJV", output_directory, ch_flag), volume, ol, 20)
        self.addToPrintList("EFC", vol(volume, volume_masks, sample_size, "EFC", output_directory, ch_flag), volume, ol, 21)
        self.addToPrintList("FBER", vol(volume, volume_masks, sample_size, "FBER", output_directory, ch_flag), volume, ol, 22)
        
    def addToPrintList(self, metric_name, value, volume, ol, il):
        
        # Add a new key-value pair to the dictionary
        self[metric_name] = value
        self["output"].append(metric_name)
        
        # Display information about the image's metrics
        if metric_name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the image with the name of <%s> is %s' % (ol, il, metric_name, volume[1], value))




#################### Image processing #######################
#############################################################


def vol(volume, volume_masks, sample_size, metric_name, output_directory, ch_flag):
    
    # Dictionary mapping each metric's name to its corresponding function
    switcher = {
            'Mean': mean,
            'Range': rang,
            'Variance': variance, 
            'CV': cv,
            'CPP': cpp,
            'PSNR': fpsnr,
            'SNR1': snr1,
            'SNR2': snr2,
            'SNR3': snr3,
            'SNR4': snr4,
            'CNR': cnr,
            'CVP': cvp,
            'CJV': cjv,
            'EFC': efc,
            'FBER': fber,
            }
    
    # Retrieve the appropriate function based on the provided metric name
    func = switcher.get(metric_name)
    Measures = []
    
    '''
    # Debug: print out volume and mask dimensions
    print(f"Volume shape: {np.shape(volume[0])}, Mask shape: {np.shape(volume_masks[0])}")
    '''
    
    # Iterate through the volume data
    for image_number, mask_number in zip(range(1, len(volume[0]), sample_size), range(1, len(volume_masks[0]), sample_size)):
        img = volume[0][image_number]
        msk = volume_masks[0][mask_number]
        
        '''
        # Debug: print out the image and mask slice indices
        print(f"Processing image index: {image_number}, mask index: {mask_number}")
        print(f"Image shape: {img.shape}, Mask shape: {msk.shape}")
        '''
        
        # Calculate foreground and background intensities
        F, B, m, f, b = foreground(img, msk, output_directory, volume, image_number)
        
        '''
        # Debug
        # Plot img, msk, fore_image, and back_image
        fig, axs = plt.subplots(1, 4, figsize=(16, 4))
        axs[0].imshow(img, cmap='gray')
        axs[0].set_title('Image (img)')
        axs[1].imshow(msk, cmap='gray')
        axs[1].set_title('Mask (msk)')
        axs[2].imshow(F, cmap='gray')
        axs[2].set_title('Foreground (fore_image)')
        axs[3].imshow(B, cmap='gray')
        axs[3].set_title('Background (back_image)')
        plt.show()
        print(f"f: {f}")
        print(f"b: {b}")
        '''
        
        # Check if foreground has zero standard deviation, skip calculation if so
        if np.std(F) == 0:
            '''
            # Debug
            print(f"Skipping slice {image_number}, no variation in foreground.")
            '''
            continue
        
        # Calculate the measure using the corresponding function
        measure = func(F, B, m, f, b)
        
        # If the measure is NaN or infinite, skip and continue
        if np.isnan(measure) or np.isinf(measure):
            continue
        
        # Append the calculated measure
        Measures.append(measure)
        
    # Return the mean of all calculated measures
    return np.mean(Measures)

    
def foreground(img, msk, output_directory, volume, image_number):
    
    # Get the image and the mask to compute the foreground and the background
    fore_image = msk * img
    back_image = (1 - msk) * img

    return fore_image, back_image, msk, img[msk == 1], img[msk == 0]


#################### Image quality control metrics ####################
#######################################################################

# All the computed measurements are average values over the entire volume, which are calculated for every slice separately.


# Mean of the foreground intensity values
def mean(F, B, m, f, b):
    return np.nanmean(f)


# Range of the foreground intensity values
def rang(F, B, m, f, b):
    if len(f) > 0:
        return np.ptp(f)
    else:
        return np.nan


# Variance of the foreground intensity values
def variance(F, B, m, f, b):
    return np.nanvar(f)


# Percentage of variation coefficent : standard deviation over the mean of the foreground intensity values
def cv(F, B, m, f, b):
    return (np.nanstd(f)/np.nanmean(f))*100


# Contrast per pixel : mean of the foreground intensity values after applying a 3x3 2D Laplacian filter
'''
def cpp(F, B, m, f, b):
    filt = np.array([[ -1/8, -1/8, -1/8],[-1/8, 1, -1/8],[ -1/8, -1/8,  -1/8]])
    I_hat = conv2(F, filt, mode='same')
    return np.nanmean(I_hat)
'''
def cpp(F, B, m, f, b):
    filt = np.array([[-1/8, -1/8, -1/8], [-1/8, 1, -1/8], [-1/8, -1/8, -1/8]])
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    I_hat = conv2(F, filt, mode='same')
    return np.nanmean(I_hat)


# Peak signal to noise ratio of the foreground intensity values
def psnr(img1, img2):
    mse = np.square(np.subtract(img1, img2)).mean()
    return 20 * np.log10(np.nanmax(img1) / np.sqrt(mse))
'''
def fpsnr(F, B, m, f, b):
    I_hat = median(F/np.max(F), square(5))
    return psnr(F, I_hat)
'''
def fpsnr(F, B, m, f, b):
    if F.ndim > 2:  # Flatten 3D arrays into 2D
        F = np.mean(F, axis=-1)
    I_hat = median(F / np.max(F), square(5))
    return psnr(F, I_hat)


# Signal to noise ratio : foreground standard deviation divided by background standard deviation
def snr1(F, B, m, f, b):
    return np.nanstd(f) / np.nanstd(b)


# Patch : random 5x5 square patch of the image
def patch(img, patch_size):
    h = int(np.floor(patch_size / 2))
    U = np.pad(img, pad_width=5, mode='constant')
    [a,b]  = np.where(img == np.max(img))
    a = a[0]
    b = b[0]
    return U[a:a+2*h+1,b:b+2*h+1]


# Signal to noise ratio : mean of the foreground patch divided by background standard deviation
def snr2(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    fore_patch = patch(F, 5)
    return np.nanmean(fore_patch) / np.nanstd(b)


# Signal to noise ratio : foreground patch standard deviation divided by the centered foreground patch standard deviation
def snr3(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    fore_patch = patch(F, 5)
    return np.nanmean(fore_patch)/np.nanstd(fore_patch - np.nanmean(fore_patch))


# Signal to noise ratio : mean of the foreground patch divided by the mean of the background patch 
def snr4(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    if B.ndim > 2:  # Flatten to 2D if F has more dimensions
        B = np.mean(B, axis=-1)
    fore_patch = patch(F, 5)
    back_patch = patch(B, 5)
    return np.nanmean(fore_patch) / np.nanstd(back_patch)


# Contrast to noise ratio : mean of the foreground and background patches difference divided by background patch standard deviation
def cnr(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    if B.ndim > 2:  # Flatten to 2D if F has more dimensions
        B = np.mean(B, axis=-1)
    fore_patch = patch(F, 5)
    back_patch = patch(B, 5)
    return np.nanmean(fore_patch-back_patch) / np.nanstd(back_patch)


# Coefficient of variation of the foreground patch : foreground patch standard deviation divided by foreground patch mean
def cvp(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    fore_patch = patch(F, 5)
    return np.nanstd(fore_patch) / np.nanmean(fore_patch)


# Coefficient of joint variation between the foreground and background
def cjv(F, B, m, f, b):
    return (np.nanstd(f) + np.nanstd(b)) / abs(np.nanmean(f) - np.nanmean(b))


# Entropy focus criterion
def efc(F, B, m, f, b):
    if F.ndim > 2:  # Flatten to 2D if F has more dimensions
        F = np.mean(F, axis=-1)
    n_vox = F.shape[0] * F.shape[1]
    efc_max = 1.0 * n_vox * (1.0 / np.sqrt(n_vox)) * \
        np.log(1.0 / np.sqrt(n_vox))
    cc = (F**2).sum()
    b_max = np.sqrt(abs(cc))
    return float((1.0 / abs(efc_max)) * np.sum(
        (F / b_max) * np.log((F + 1e16) / b_max)))


# Foreground-background energy ratio
def fber(F, B, m, f, b):
    fg_mu = np.nanmedian(np.abs(f) ** 2)
    bg_mu = np.nanmedian(np.abs(b) ** 2)
    if bg_mu < 1.0e-3:
        return 0
    return float(fg_mu / bg_mu)


# Slice orientation
def orientation(scan):
    image_data, image_header = load(scan)
    image_shape = np.shape(image_data)
    if image_shape[0] < image_shape[1] and image_shape[0] < image_shape[2]:
        orien = 0
    elif image_shape[1] < image_shape[0] and image_shape[1] < image_shape[2]:
        orien = 1
    elif image_shape[2] < image_shape[0] and image_shape[2] < image_shape[1]:
        orien = 2
    else:
        orien = 3
    return orien


#################### CSV file editing ####################
##########################################################


def worker_callback(measures, output_directory):
    
    # Access global variables
    global csv_report, first, nfiledone
    
    # Check that the first file is being processed
    if nfiledone  == 0:
        # Open the CSV report file in append mode or overwrite mode
        csv_report = open(output_directory + os.sep + "results" + ".tsv" , overwrite_flag, buffering=1)
        first = True

    if first and overwrite_flag == "w": 
        first = False
        # Write comment lines for headers
        csv_report.write("\n".join(["#" + measures for measures in headers])+"\n")
        csv_report.write("#dataset:"+"\t".join(measures["output"])+"\n")
    
    # Write data to the CSV report file and replaces nan values by 0.0
    # csv_report.write("\t".join([str(s[field]) for field in s["output"]])+"\n")          # replace replace_nan by str
    csv_report.write("\t".join([str(0.0) if isinstance(measures[field], float) and math.isnan(measures[field]) else str(measures[field]) for field in measures["output"]]) + "\n")  

    csv_report.flush()      # Flush the buffer to ensure writing immediately
    nfiledone += 1          # Increment the count of processed files
    print('The results are updated.')


#################### Data dimensionality reduction and classification ####################
##########################################################################################


def tsv_to_dataframe(tsvfileaddress):       # Creates a dataframe with all the extracted metadata and calculated data from images
    # This creates the dataframe used for the tables, graphs, charts, UMAP and TSNE.
    # If you change things here (columns or rows used), it will modify what is shown in the user interface
    # Read the CSV file into a pandas dataframe, skipping the first two rows and using the third row as headers
    return pd.read_csv(tsvfileaddress, sep='\t', skiprows=2, header=0)


def data_whitening(dframe):
    
    # Fill missing values with N/A in the dataframe
    dframe = dframe.fillna('N/A')
    
    # Create a copy of the dataframe     
    df = dframe.copy()
    
    # Select columns excluding object type
    df = df.select_dtypes(exclude=['object'])
    
    # Enter the names of the columns you want included for the UMAP and TSNE 
    df = df[["MEAN", "RNG", "VAR", "CV", "CPP", "PSNR","SNR1", "SNR2", "SNR3", "SNR4", "CNR", "CVP", "CJV", "EFC", "FBER"]]
    
    # Apply whitening transformation to the selected numeric data. The data are transformed to be used by UMAP and TSNE
    ds = whiten(df)  
    
    # Returns the transformed dataset : these are the data used by UMAP and TSNE
    return ds


def tsne_umap(dataframe, per):
    
    # Apply data whitening to the dataframe. Is used for TSNE and UMAP
    ds = data_whitening(dataframe)
    
    # Performs t-SNE on the transformed dataset. Perplexity must be less than n subjects. default = 30. Usually between 5 and 50
    tsne = TSNE(n_components=2, random_state=0, perplexity = 5)
    
    # Learns the parameters, fit to data and applies the transformation to new dataset
    tsne_obj = tsne.fit_transform(ds)       
    
    # Add t-SNE components to the original DataFrame as 'x' and 'y'
    dataframe['x'] = tsne_obj[:,0].astype(float)
    dataframe['y'] = tsne_obj[:,1].astype(float)
    
    # Performs UMAP dimensionality reduction
    reducer = umap.UMAP()
    
    # Learns the parameters, fit to data and applies the transformation to new dataset                               
    reducer_obj = reducer.fit_transform(ds)     
    
    # Add UMAP components to the original DataFrame as 'u' and 'v'
    dataframe['u'] = reducer_obj[:,0]
    dataframe['v'] = reducer_obj[:,1]


def cleanup(final_address, per):
    
    # Read data from the final TSV file into a dataframe
    df = tsv_to_dataframe(final_address)
    
    # Apply t-SNE and UMAP dimensionality reduction on the dataframe
    tsne_umap(df, per)
    
    # Read only the header row from the final TSV file
    hf = pd.read_csv(final_address, sep='\t',  nrows=1)
    
    # Write the header row back to the same CSV file (overwriting the existing file)
    hf.to_csv(final_address, index = None, header=True, sep = '\t', mode = 'w')
    
    # Append the modified dataframe to the CSV file (below the header)
    df.to_csv(final_address, index = None, header=True, sep = '\t', mode = 'a')
    
    # Return the modified dataframe
    return df


#############################################################################################################################
##########################              Main function calling all the other functions              ##########################
#############################################################################################################################


if __name__ == '__main__':
    
    # Starting message
    print("MRQy is starting ...")
    
    # Record the start time for runtime measurement
    start_time = time.time()
    
    # Add the start time information to the headers list
    headers.append(f"start_time:\t{datetime.datetime.now()}")
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description = '')
    
    # Path for the output folder name
    parser.add_argument('output_folder_name',
                        help = "the subfolder name on the '...\\UserInterface\\Data\\output_folder_name' directory.",
                        type=str)
    
    # Path for the input directory
    parser.add_argument('inputdir',
                        help = "input foldername consists of *.nii.gz files. For example: 'E:\\Data\\Rectal\\input_data_folder'",
                        nargs = "*")
    
    # Other command line arguments
    parser.add_argument('-r', help="folders as name", default=False)
    parser.add_argument('-s', help="save foreground masks", default=False)
    parser.add_argument('-b', help="number of samples", default=1, type = int)
    parser.add_argument('-u', help="percent of middle images", default=100)
    parser.add_argument('-t', help="the address of the user-specified tags list (*.txt)", default=0)
    parser.add_argument('-c', help="if yes the ch computes objects", default=False)
    
    args = parser.parse_args() 
    root = args.inputdir[0]         # root = input directory
    
    # Setting flags based on parsed arguments
    if args.r == 0:
        folders_flag = "False"
    else: 
        folders_flag = args.r
    if args.s == 0:
        save_masks_flag = "False" 
    else: 
        save_masks_flag = args.s
    if args.b != 1:
        sample_size = args.b
    else: 
        sample_size = 1
    if args.u != 100:
        middle_size = int(args.u)
    else: 
        middle_size = 100
    if args.t != 0:
        tag_names = [line.strip() for line in open(args.t, "r")]
        tag_list = [line.strip().replace(" ", "") for line in open(args.t, "r")] 
    if args.c == 0:
        ch_flag = "False"
    else: 
        ch_flag = args.c
    
    # The results are saved in the "UserInterface" directory
    print_folder_note = os.getcwd() + os.sep + 'UserInterface'
    
    # A "Data" folder in UserInterface will contain the analysis results in the subfolder "output_folder_name"
    # output_directory = UserInterface/Data/output_folder_name
    output_directory = print_folder_note + os.sep + 'Data' + os.sep + args.output_folder_name
    
    overwrite_flag = "w"        
    headers.append(f"outdir:\t{os.path.realpath(output_directory)}")
    
    # Function file_name(root) : all the files 
    files, nifti_names, nifti_images = file_name(root)                                               
    
    # Determine data types and set corresponding flags (if there are niftis or not)
    if len(nifti_images) > 0:
        nifti_flag = True
    if len(nifti_images) == 0:
        nifti_flag = False
        print('The input folder is empty or includes unsupported files format!')
    
    output_masks = os.path.join(output_directory, 'masks')
    os.makedirs(output_masks, exist_ok=True)

    # Process NIFTI data if available
    for i in range(len(nifti_names)):
        if nifti_flag:
            for l, scan in enumerate(nifti_images):
                mask = brain_extraction(scan, output_masks)
                volume = volume_nifti(scan, nifti_names[l])
                volume_masks = volume_nifti_masks(mask)
                saveThumbnails_nifti(volume, output_directory)
                # saveThumbnails_nifti_masks(volume_masks, output_directory)
                measures = BaseVolume_nifti(output_directory, volume, l+1, scan, sample_size, ch_flag)
                # contours = contours(scan, mask)
                worker_callback(measures, output_directory)
            nifti_flag = False
    
    # Create the path for the result TSV file
    address = output_directory + os.sep + "results" + ".tsv" 
    
    # Process data for dimensionality reduction
    if len(nifti_names) < 6:
        # t-SNE and UMPA cannot be performed if we have less than 6 images
        print('Skipped the t-SNE and UMAP computation because of insufficient data. The UMAP and t-SNE process need at least 6 input data.')
        df = tsv_to_dataframe(address)
    
    else:    
        df = cleanup(address, 30)
        df = df.drop(['Name of Images'], axis=1)
        df = df.rename(columns={"#dataset:Image": "Image", "x":"TSNEX","y":"TSNEY", "u":"UMAPX", "v":"UMAPY" })
        df = df.fillna('N/A')
    
    # Save processed data as IQM.csv
    df.to_csv(output_directory + os.sep + 'IQM.csv', index=False)
    print("The IQMs data are saved in the {} file. ".format(output_directory + os.sep + "IQM.csv"))
    
    # Print execution information
    print("Done!")
    print("MRQy program took", format((time.time() - start_time)/60, '.2f'), \
          "minutes for {} images to run.".format(len(nifti_names)))
    
    # Provide guidance for viewing the final results in the MRQy interface
    msg = "Please go to the '{}' directory and open up the 'index.html' file.\n".format(print_folder_note) + \
    "Click on 'View Results' and select '{}' file.\n".format(output_directory + os.sep + "results.tsv")   
    print_msg_box(msg, indent=3, width=None, title="To view the final MRQy interface results:")


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################