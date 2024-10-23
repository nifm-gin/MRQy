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
from itertools import accumulate                            # Return series of accumulated sums
import matplotlib.cm as cm                                  # Interactive colormaps
import matplotlib.pyplot as plt                             # Interactive plots
from medpy.io import load                                   # For NIFTI image processing
import numpy as np                                          # Mathematical operations over arrays
import os                                                   # Interaction with the operating system
import pandas as pd                                         # Data analysis
import pydicom                                              # Read DICOM files
from scipy.io import loadmat                                # Load MATLAB file
from scipy.ndimage import rotate                            # Rotation of arrays
from scipy.signal import convolve2d as conv2                # Signal processing, convolution
from skimage.filters import median                          # Image local median
from skimage.measure import find_contours                   # Find iso-valued contours in a 2D array for a given level value
from skimage.morphology import square                       # Square shaped footprint of an image
import subprocess                                           # use command-line in python
import time                                                 # Manipulate time values
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

#############################################################################################################################
##############################                           Definitions                           ##############################
#############################################################################################################################


#################### File pre-processing ####################
#############################################################


def file_name(root):
    
    # Starting message
    print('MRQy is starting....')
    
    # Gathering relevant file paths based on extensions. MRQy supports .dcm, .nii.gz and .mat files
    files = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(root) 
                for filename in filenames 
                if filename.endswith('.dcm') 
                or filename.endswith('.nii.gz')
                or filename.endswith('.mat')]
    
    # Separating files based on their extensions
    mats = [i for i in files if i.endswith('.mat')]
    dicoms = [i for i in files if i.endswith('.dcm')]
    niftis = [i for i in files if i.endswith('.nii.gz')]

    
    # Extracting subject identifiers from the different files
    niftis_subjects = [os.path.basename(scan)[:os.path.basename(scan).index('.')] for scan in niftis]
    dicom_subjects = []
    mat_subjects = [os.path.basename(scan)[:os.path.basename(scan).index('.')] for scan in mats]
    
    # For DICOMS only
    if folders_flag == "False":
        # Processing individual DICOM files
        for i in dicoms: 
            dicom_subjects.append(pydicom.dcmread(i).PatientID) 
        # Counting occurences of each patient ID
        duplicateFrequencies = {}
        for i in dicom_subjects:
            duplicateFrequencies[i] = dicom_subjects.count(i)
        # Extracting unique patient IDs and their counts
        subjects_id = []
        subjects_number = []
        for i in range(len(duplicateFrequencies)):
              subjects_id.append(list(duplicateFrequencies.items())[i][0])
              subjects_number.append(list(duplicateFrequencies.items())[i][1])
        ind = [0] + list(accumulate(subjects_number))
        # Splitting DICOM files within folders
        splits = [dicoms[ind[i]:ind[i+1]] for i in range(len(ind)-1)]
    elif folders_flag == "True":
        # Processing DICOM files within folders
        dicom_subjects = [d for d in os.listdir(root) if os.path.isdir(root + os.sep + d)]
        subjects_number = []
        for i in range(len(dicom_subjects)):
            # Counting DICOM files per subject
            subjects_number.append(
                len([os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(root + os.sep + dicom_subjects[i]) 
            for filename in filenames 
            if filename.endswith('.dcm')]))
        subjects_id  = dicom_subjects
        ind = [0] + list(accumulate(subjects_number))
        # Splitting DICOM files based on folders (subjects)
        splits = [dicoms[ind[i]:ind[i+1]] for i in range(len(ind)-1)]

    # Combining all subject identifiers
    subjects = subjects_id + niftis_subjects + mat_subjects
    # Displaying the total number of identified subjects
    print('The number of images is {}'.format(len(subjects)))
    # Returning various lists containing file paths, subjects, and DICOM splits
    return files, subjects, splits, niftis, niftis_subjects, mats, mat_subjects


#################### DICOM files processing ####################
################################################################


def volume_dicom(scans, name):
    # Selecting a portion of scans based on size
    scans = scans[int(0.005 *len(scans)*(100 - middle_size)):int(0.005 *len(scans)*(100 + middle_size))]
    # Reading metadata from the first DICOM file
    inf = pydicom.dcmread(scans[0])
    # Modifying attributes if they exist         
    if hasattr(inf, 'MagneticFieldStrength'):
        if inf.MagneticFieldStrength > 10:
            inf.MagneticFieldStrength = inf.MagneticFieldStrength/10000
    else:
        inf.MagneticFieldStrength = ''
    if hasattr(inf, 'Manufacturer') == False:
        inf.Manufacturer = ''
    if  hasattr(inf, 'RepetitionTime') == False:
            inf.RepetitionTime = 0
    if  hasattr(inf, 'EchoTime') == False:
            inf.EchoTime = 0
    # Determining name value based on folders_flag
    if folders_flag == "False":
        name_value = inf.PatientID
    elif folders_flag == "True":
        name_value = name
    # Slice orientation
    orientation = inf.ImageOrientationPatient
    row_cosines = orientation[:3]
    col_cosines = orientation[3:]
    if abs(row_cosines[0]) > 0.9:
        ori = 0         # Sagittal
        return ori    
    elif abs(row_cosines[1]) > 0.9:
        ori = 1         # Coronal
        return ori    
    elif abs(row_cosines[2]) > 0.9:
        ori = 2         # Axial
        return ori             
    # Creating a dictionnary of DICOM metadata attributes and their values    
    tags = {
             'ID': name_value,                                  # Image / Subject ID
             'MFR': inf.Manufacturer,                           # Manufacturer name from the file header
             'VRX': format(inf.PixelSpacing[0], '.2f'),         # Voxel resolution in x plane
             'VRY': format(inf.PixelSpacing[1], '.2f'),         # Voxel resolution in y plane
             'VRZ': format(inf.SliceThickness, '.2f'),          # Voxel resolution in z plane
             'MFS': inf.MagneticFieldStrength,                  # Magnetic fiels strength from the file header
             'ROWS': int(inf.Rows),                             # Rows value of the volume
             'COLS': int(inf.Columns),                          # Columns value of the volume
             'TR': format(inf.RepetitionTime, '.2f'),           # Repetition time value of the volume
             'TE': format(inf.EchoTime, '.2f'),                 # Echo time value of the volume
             'NUM': len(scans),                                 # Number of slice images in each volume
             'ORIENTATION': ori                                 # Slice orientation : Sagittal (0), Coronal (1) or Axial (2)
    }
    # Fetching additional attributes if available
    tag_values = []
    # Checking for args variable 
    if args.t != 0:
        for de in tag_list:
            # if hasattr(inf, de) == False or inf.data_element(de).value == '':
            if hasattr(inf, de) == False:
                value = ''
            else:
                value = inf.data_element(de).value
            tag_values.append(value)
        res_dct = dict(zip(iter(tag_names), iter(tag_values)))
        tags.update(res_dct)
    # Reading and sorting DICOM files, creating a 3D image volume    
    slices = [pydicom.read_file(s) for s in scans]
    slices.sort(key = lambda x: int(x.InstanceNumber))
    # PL = pd.DataFrame([s.pixel_array for s in slices], columns=['images'])
    # images = PL['images'].to_numpy().astype(np.int64)
    images = np.stack([s.pixel_array for s in slices])
    images = images.astype(np.int64)
    # Returning the image volume and associated tags
    return images, tags


def saveThumbnails_dicom(v, output):
    # Check if saving masks is enabled
    if save_masks_flag!='False':
        ffolder = output + '_foreground_masks'
        # Create a directory for foreground masks
        os.makedirs(ffolder + os.sep + v[1]['ID'])
    elif save_masks_flag=='False':
        ffolder = output
    # Create a directory for images
    # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    os.makedirs(output + os.sep + v[1]['ID'])
    # Save images as thumbnails
    for i in range(0, len(v[0]), sample_size):
        # Le chemin + nom qu'ont les images.png dans leur dossier respectif
        plt.imsave(output + os.sep + v[1]['ID'] + os.sep + v[1]['ID'] + '(%d).png' % i, v[0][i], cmap = cm.Greys_r)
    # Print the number of saved images and the directory
    print('The number of %d images are saved to %s' % (len(v[0]),output + os.sep + v[1]['ID']))
    return ffolder + os.sep + v[1]['ID']


class BaseVolume_dicom(dict):

    def __init__(self, fname_outdir, v, ol, folder_foregrounds, sample_size, ch_flag):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        # Add image information to the output list
        self.addToPrintList("Image", v[1]['ID'], v, ol, 170)
        # Add the directory content to the output list
        self["outdir"] = fname_outdir
        self.addToPrintList("Name of Images", os.listdir(fname_outdir + os.sep + v[1]['ID']), v, ol, 100)
        # Add specific data attributes to the output list
        for i,j in enumerate(v[1]):
            if i != 0:
                self.addToPrintList(j, v[1][j], v, ol, i)
        
        # Add information to the ouptut list
        self.addToPrintList("MFR", v[1]['Manufacturer'], v, ol, 1)
        self.addToPrintList("MFS", v[1]['MFS'], v, ol, 2)
        self.addToPrintList("VRX", v[1]['VR_x'], v, ol, 3)
        self.addToPrintList("VRY", v[1]['VR_y'], v, ol, 4)
        self.addToPrintList("VRZ", v[1]['VR_z'], v, ol, 5)
        self.addToPrintList("ROWS", v[1]['Rows'], v, ol, 6)
        self.addToPrintList("COLS", v[1]['Columns'], v, ol, 7)
        self.addToPrintList("TR", v[1]['TR'], v, ol, 8)
        self.addToPrintList("TE", v[1]['TE'], v, ol, 9)
        self["os_handle"] = v[0]
        self.addToPrintList("NUM", v[1]['Number'], v, ol, 10)
        self.addToPrintList("ORIENTATION", v[1]['Orientation'], v, ol, 11)
        self.addToPrintList("MEAN", vol(v, sample_size, "Mean", folder_foregrounds, ch_flag), v, ol, 12)
        self.addToPrintList("RNG", vol(v, sample_size, "Range", folder_foregrounds, ch_flag), v, ol, 13)
        self.addToPrintList("VAR", vol(v, sample_size, "Variance", folder_foregrounds, ch_flag), v, ol, 14)
        self.addToPrintList("CV", vol(v, sample_size, "CV", folder_foregrounds, ch_flag), v, ol, 15)
        self.addToPrintList("CPP", vol(v, sample_size, "CPP", folder_foregrounds, ch_flag), v, ol, 16)
        self.addToPrintList("PSNR", vol(v, sample_size, "PSNR", folder_foregrounds, ch_flag), v, ol, 17)
        self.addToPrintList("SNR1", vol(v, sample_size, "SNR1", folder_foregrounds, ch_flag), v, ol, 18)
        self.addToPrintList("SNR2", vol(v, sample_size, "SNR2", folder_foregrounds, ch_flag), v, ol, 19)
        self.addToPrintList("SNR3", vol(v, sample_size, "SNR3", folder_foregrounds, ch_flag), v, ol, 20)
        self.addToPrintList("SNR4", vol(v, sample_size, "SNR4", folder_foregrounds, ch_flag), v, ol, 21)
        self.addToPrintList("CNR", vol(v, sample_size, "CNR", folder_foregrounds, ch_flag), v, ol, 22)
        self.addToPrintList("CVP", vol(v, sample_size, "CVP", folder_foregrounds, ch_flag), v, ol, 23)
        self.addToPrintList("CJV", vol(v, sample_size, "CJV", folder_foregrounds, ch_flag), v, ol, 24)
        self.addToPrintList("EFC", vol(v, sample_size, "EFC", folder_foregrounds, ch_flag), v, ol, 25)
        self.addToPrintList("FBER", vol(v, sample_size, "FBER", folder_foregrounds, ch_flag), v, ol, 26)
        
    def addToPrintList(self, name, val, v, ol, il):
        # Add a new key-value pair to the dictionary
        self[name] = val
        self["output"].append(name)
        # Display information about the image's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the image with the name of <%s> is %s' % (ol,il,name, v[1]['ID'], val))


#################### NIFTI files processing ####################
################################################################


def brain_extraction(scan, output_masks):
    # Make the path for the output masks
    output_mask = os.path.join(output_masks, os.path.basename(scan))
    
    # FSL-BET
    # subprocess.run(['bet', scan, output_mask, '-f', '0.5', '-g', '0', '-m', '-n'], check=True)
    
    # HD-BET fast 
    subprocess.run(['hd-bet', '-i', scan, '-o', output_mask, '-tta', '0', '-pp', '0', '-b', '0', '-mode', 'fast'], check=True)
    
    # HD-BET accurate
    # subprocess.run(['hd-bet', '-i', scan, '-o', output_mask, '-tta', '1', '-pp', '1', '-b', '0', '-mode', 'accurate'], check=True)
    
    mask_path = output_mask.replace(".nii.gz", "_mask.nii.gz")

    return mask_path


def volume_nifti(scan, name):
    
    # Loading image data and header
    image_data, image_header = load(scan)
    # Get image dimensions
    image_shape = np.shape(image_data)
    # Extract the 2D images from the 3D image data
    if image_shape[0] < image_shape[1] and image_shape[0] < image_shape[2]:         # Sagittal
        images = [image_data[i, :, :] for i in range(image_shape[0])]
    elif image_shape[1] < image_shape[0] and image_shape[1] < image_shape[2]:       # Coronal
        images = [image_data[:, i, :] for i in range(image_shape[1])]               
    elif image_shape[2] < image_shape[0] and image_shape[2] < image_shape[1]:       # Axial
        images = [image_data[:, :, i] for i in range(image_shape[2])]
    else:                                                                           # 3D
        images = [image_data[:, :, i] for i in range(image_shape[2])]
    
    # Return image, name and image header
    return images, name, image_header


def volume_nifti_masks(mask_path):
    
    # Transform mask data and header
    mask_data, mask_header = load(mask_path)
    # Get mask dimensions
    mask_shape = np.shape(mask_data)
    # Extract the 2D masks from the 3D mask data
    if mask_shape[0] < mask_shape[1] and mask_shape[0] < mask_shape[2]:         # Sagittal
        masks = [mask_data[i, :, :] for i in range(mask_shape[0])]
    elif mask_shape[1] < mask_shape[0] and mask_shape[1] < mask_shape[2]:       # Coronal
        masks = [mask_data[:, i, :] for i in range(mask_shape[1])]               
    elif mask_shape[2] < mask_shape[0] and mask_shape[2] < mask_shape[1]:       # Axial
        masks = [mask_data[:, :, i] for i in range(mask_shape[2])]
    else:                                                                           # 3D
        masks = [mask_data[:, :, i] for i in range(mask_shape[2])]
        
    mask_name = os.path.basename(mask_path).split('.')[0]
    
    # Return image, name and image header
    return masks, mask_name, mask_header


def saveThumbnails_nondicom(v, output, masks):
    # Create a directory for images
    # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    os.makedirs(output + os.sep + v[1])
    # Save images as thumbnails, with the mask contours
    for i in range(len(v[0])):
        # Extract the image and the corresponding mask slice
        img_slice = v[0][i]
        mask_slice = masks[i]
        # Find contours from the mask
        img_slice = rotate(img_slice, 90)
        mask_slice = rotate(mask_slice, 90)
        contours = find_contours(mask_slice, level=0.5)
        # Plot image and overlay contours
        fig, ax = plt.subplots()
        ax.imshow(img_slice, cmap='gray')
        # Plot the contours of the ROI
        for contour in contours:
            ax.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
        # Plot the contours of the images with a mask
        if np.any(mask_slice != 0):
            ax.plot([0, img_slice.shape[1], img_slice.shape[1], 0, 0], [0, 0, img_slice.shape[0], img_slice.shape[0], 0],
                    color='green', linewidth=10)
        ax.axis('off')  # Turn off the axis
        plt.savefig(output + os.sep + v[1] + os.sep + v[1] + f'({i+1}).png', bbox_inches='tight', pad_inches=0)
        plt.close(fig)  # Close the plot to save memory    
    # print('image number %d out of %d is saved to %s' % (int(i+1), len(v[0]),output + os.sep + v[1]))
    print('The number of %d images are saved to %s' % (len(v[0]),output + os.sep + v[1]))


class BaseVolume_nondicom(dict):

    def __init__(self, fname_outdir, v, ol, scan, sample_size, ch_flag):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        
        # Add information to the output list
        self.addToPrintList("Image", v[1], v, ol, 170)
        self["outdir"] = fname_outdir
        self.addToPrintList("Name of Images", os.listdir(fname_outdir + os.sep + v[1]), v, ol, 100)
        self.addToPrintList("VRX", format(v[2].get_voxel_spacing()[0], '.2f'), v, ol, 1)
        self.addToPrintList("VRY", format(v[2].get_voxel_spacing()[1], '.2f'), v, ol, 2)
        self.addToPrintList("VRZ", format(v[2].get_voxel_spacing()[2], '.2f'), v, ol, 3)
        self.addToPrintList("ROWS", np.shape(v[0])[1], v, ol, 4)
        self.addToPrintList("COLS", np.shape(v[0])[2], v, ol, 5)
        self["os_handle"] = v[0]
        self.addToPrintList("NUM", len(v[0]), v, ol, 6)
        self.addToPrintList("ORIENTATION", orientation(scan), v, ol, 7) 
        self.addToPrintList("MEAN", vol(v, volume_masks, sample_size, "Mean", fname_outdir, ch_flag), v, ol, 8)
        self.addToPrintList("RNG", vol(v, volume_masks, sample_size, "Range", fname_outdir, ch_flag), v, ol, 9)
        self.addToPrintList("VAR", vol(v, volume_masks, sample_size, "Variance", fname_outdir, ch_flag), v, ol, 10)
        self.addToPrintList("CV", vol(v, volume_masks, sample_size, "CV", fname_outdir, ch_flag), v, ol, 11)
        self.addToPrintList("CPP", vol(v, volume_masks, sample_size, "CPP", fname_outdir, ch_flag), v, ol, 12)
        self.addToPrintList("PSNR", vol(v, volume_masks, sample_size, "PSNR", fname_outdir, ch_flag), v, ol, 13)
        self.addToPrintList("SNR1", vol(v, volume_masks, sample_size, "SNR1", fname_outdir, ch_flag), v, ol, 14)
        self.addToPrintList("SNR2", vol(v, volume_masks, sample_size, "SNR2", fname_outdir, ch_flag), v, ol, 15)
        self.addToPrintList("SNR3", vol(v, volume_masks, sample_size, "SNR3", fname_outdir, ch_flag), v, ol, 16)
        self.addToPrintList("SNR4", vol(v, volume_masks, sample_size, "SNR4", fname_outdir, ch_flag), v, ol, 17)
        self.addToPrintList("CNR", vol(v, volume_masks, sample_size, "CNR", fname_outdir, ch_flag), v, ol, 18)
        self.addToPrintList("CVP", vol(v, volume_masks, sample_size, "CVP", fname_outdir, ch_flag), v, ol, 19)
        self.addToPrintList("CJV", vol(v, volume_masks, sample_size, "CJV", fname_outdir, ch_flag), v, ol, 20)
        self.addToPrintList("EFC", vol(v, volume_masks, sample_size, "EFC", fname_outdir, ch_flag), v, ol, 21)
        self.addToPrintList("FBER", vol(v, volume_masks, sample_size, "FBER", fname_outdir, ch_flag), v, ol, 22)
        
    def addToPrintList(self, name, val, v, ol, il):
        # Add a new key-value pair to the dictionary
        self[name] = val
        self["output"].append(name)
        # Display information about the image's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the image with the name of <%s> is %s' % (ol,il,name, v[1], val))


##################### MAT files processing #####################
################################################################


def volume_mat(mat_scan, name):
    # Loading volume data 
    v1 = loadmat(mat_scan)['vol']
    # Creating a dictionnary for the ID
    tags = {'ID': name}
    # Return the loaded volume data and the tags
    return v1, tags


def saveThumbnails_mat(v, output):
    # Check if saving masks is enabled
    if save_masks_flag!='False':
        ffolder = output + '_foreground_masks'
        # Create a directory for foreground masks
        os.makedirs(ffolder + os.sep + v[1]['ID'])
    elif save_masks_flag=='False':
        ffolder = output
    # Create a directory for images 
    # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    os.makedirs(output + os.sep + v[1]['ID'])
    # Save image as thumbnails
    for i in range(np.shape(v[0])[2]):
        # Le chemin + nom qu'ont les images.png dans leur dossier respectif
        plt.imsave(output + os.sep + v[1]['ID']+ os.sep + v[1]['ID'] + '(%d).png' % int(i+1), v[0][:,:,i], cmap = cm.Greys_r)
        # Print the number of saved images and the directory
    print('The number of %d images are saved to %s' % (np.shape(v[0])[2],output + os.sep + v[1]['ID']))
    return ffolder + os.sep + v[1]['ID']


class BaseVolume_mat(dict):

    def __init__(self, fname_outdir, v, ol,folder_foregrounds, sample_size):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        
        # Add image information to the output list
        self.addToPrintList("Image", v[1]['ID'], v, ol, 170)
        
         # Add information to the output list
        self["outdir"] = fname_outdir
        self.addToPrintList("Name of Images", os.listdir(fname_outdir + os.sep + v[1]['ID']), v, ol, 100)
        self.addToPrintList("ROWS", np.shape(v[0])[0], v, ol, 1)
        self.addToPrintList("COLS", np.shape(v[0])[1], v, ol, 2)
        self["os_handle"] = v[0]
        self.addToPrintList("NUM", np.shape(v[0])[2], v, ol, 3)
        self.addToPrintList("MEAN", vol(v, sample_size, "Mean", folder_foregrounds), v, ol, 4)
        self.addToPrintList("RNG", vol(v, sample_size, "Range", folder_foregrounds), v, ol, 5)
        self.addToPrintList("VAR", vol(v, sample_size, "Variance", folder_foregrounds), v, ol, 6)
        self.addToPrintList("CV", vol(v, sample_size, "CV", folder_foregrounds), v, ol, 7)
        self.addToPrintList("CPP", vol(v, sample_size, "CPP",folder_foregrounds), v, ol, 8)
        self.addToPrintList("PSNR", vol(v, sample_size, "PSNR", folder_foregrounds), v, ol, 9)
        self.addToPrintList("SNR1", vol(v, sample_size, "SNR1", folder_foregrounds), v, ol, 10)
        self.addToPrintList("SNR2", vol(v, sample_size, "SNR2", folder_foregrounds), v, ol, 11)
        self.addToPrintList("SNR3", vol(v, sample_size, "SNR3", folder_foregrounds), v, ol, 12)
        self.addToPrintList("SNR4", vol(v, sample_size, "SNR4", folder_foregrounds), v, ol, 13)
        self.addToPrintList("CNR", vol(v, sample_size, "CNR", folder_foregrounds), v, ol, 14)
        self.addToPrintList("CVP", vol(v, sample_size, "CVP", folder_foregrounds), v, ol, 15)
        self.addToPrintList("CJV", vol(v, sample_size, "CJV", folder_foregrounds), v, ol, 16)
        self.addToPrintList("EFC", vol(v, sample_size, "EFC", folder_foregrounds), v, ol, 17)
        self.addToPrintList("FBER", vol(v, sample_size, "FBER", folder_foregrounds), v, ol, 18)
        
    def addToPrintList(self, name, val, v, ol, il):
        # Add a new key-value pair to the dictionary
        self[name] = val
        self["output"].append(name)
        # Display information about the image's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the image with the name of <%s> is %s' % (ol,il,name, v[1]['ID'], val))


#################### Image processing #######################
#############################################################


def vol(v, volume_masks, sample_size, kk, outi_folder, ch_flag):
    # Dictionary mapping each metric's name to its corresponding function
    switcher={
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
    func=switcher.get(kk)
    M = []
    # Iterate through the volume data
    for i, mask_number in zip(range(1, len(v[0]), sample_size), range(1, len(volume_masks[0]), sample_size)):
        I = v[0][i]
        msk = volume_masks[0][mask_number]
        # If the mask is null, then ignore slice
        if np.all(msk == 0):
            continue
        # Calculate foreground and background intensities
        F, B, c, f, b = foreground(I, msk, outi_folder, v, i)
        # Check if the standard deviation of foreground is zero, skip computing measures
        if np.std(F) == 0:  # whole zero slice, no measure computing
            continue
        # Calculate the measure using the corresponding function
        measure = func(F, B, c, f, b)
        # If the measure is NaN or infinite, skip and continue
        if np.isnan(measure) or np.isinf(measure):
            continue
        # Append the calculated measure
        M.append(measure)
    # Return the mean of all calculated measures
    return np.mean(M)


def foreground(img, msk, save_folder, v, inumber):
    
    # Get the image and the mask to compute the foreground and the background
    fore_image = msk*img
    back_image = (1 - msk) * img
         
    # if not os.path.isdir(save_folder + os.sep + v[1]['ID']):                              
    return fore_image, back_image, msk, img[msk == 1], img[msk == 0]
    # fore_image = F    back_image = B      conv_hull = c       img[conv_hull] = f      img[conv_hull==False] = b


#################### Image quality control metrics ####################
#######################################################################

# All the computed measurements are average values over the entire volume, which are calculated for every slice


# Mean of the foreground intensity values
def mean(F, B, c, f, b):
    return np.nanmean(f)


# Range of the foreground intensity values
def rang(F, B, c, f, b):
    if len(f) > 0:
        return np.ptp(f)
    else:
        return np.nan


# Variance of the foreground intensity values
def variance(F, B, c, f, b):
    return np.nanvar(f)


# Percentage of variation coefficent : standard deviation over the mean of the foreground intensity values
def cv(F, B, c, f, b):
    return (np.nanstd(f)/np.nanmean(f))*100


# Contrast per pixel : mean of the foreground intensity values after applying a 3x3 2D Laplacian filter
def cpp(F, B, c, f, b):
    filt = np.array([[ -1/8, -1/8, -1/8],[-1/8, 1, -1/8],[ -1/8, -1/8,  -1/8]])
    I_hat = conv2(F, filt, mode='same')
    return np.nanmean(I_hat)


# Peak signal to noise ratio of the foreground intensity values
def psnr(img1, img2):
    mse = np.square(np.subtract(img1, img2)).mean()
    return 20 * np.log10(np.nanmax(img1) / np.sqrt(mse))
def fpsnr(F, B, c, f, b):
    I_hat = median(F/np.max(F), square(5))
    return psnr(F, I_hat)


# Signal to noise ratio : foreground standard deviation divided by background standard deviation
def snr1(F, B, c, f, b):
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
def snr2(F, B, c, f, b):
    fore_patch = patch(F, 5)
    return np.nanmean(fore_patch) / np.nanstd(b)


# Signal to noise ratio : foreground patch standard deviation divided by the centered foreground patch standard deviation
def snr3(F, B, c, f, b):
    fore_patch = patch(F, 5)
    return np.nanmean(fore_patch)/np.nanstd(fore_patch - np.nanmean(fore_patch))


# Signal to noise ratio : mean of the foreground patch divided by the mean of the background patch 
def snr4(F, B, c, f, b):
    fore_patch = patch(F, 5)
    back_patch = patch(B, 5)
    return np.nanmean(fore_patch) / np.nanstd(back_patch)


# Contrast to noise ratio : mean of the foreground and background patches difference divided by background patch standard deviation
def cnr(F, B, c, f, b):
    fore_patch = patch(F, 5)
    back_patch = patch(B, 5)
    return np.nanmean(fore_patch-back_patch) / np.nanstd(back_patch)


# Coefficient of variation of the foreground patch : foreground patch standard deviation divided by foreground patch mean
def cvp(F, B, c, f, b):
    fore_patch = patch(F, 5)
    return np.nanstd(fore_patch) / np.nanmean(fore_patch)


# Coefficient of joint variation between the foreground and background
def cjv(F, B, c, f, b):
    return (np.nanstd(f) + np.nanstd(b)) / abs(np.nanmean(f) - np.nanmean(b))


# Entropy focus criterion
def efc(F, B, c, f, b):
    n_vox = F.shape[0] * F.shape[1]
    efc_max = 1.0 * n_vox * (1.0 / np.sqrt(n_vox)) * \
        np.log(1.0 / np.sqrt(n_vox))
    cc = (F**2).sum()
    b_max = np.sqrt(abs(cc))
    return float((1.0 / abs(efc_max)) * np.sum(
        (F / b_max) * np.log((F + 1e16) / b_max)))


# Foreground-background energy ratio
def fber(F, B, c, f, b):
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


def worker_callback(s,fname_outdir):
    # Access global variables
    global csv_report, first, nfiledone
    # Check that the first file is being processed
    if nfiledone  == 0:
        # Open the CSV report file in append mode or overwrite mode
        csv_report = open(fname_outdir + os.sep + "results" + ".tsv" , overwrite_flag, buffering=1)
        first = True

    if first and overwrite_flag == "w": 
        first = False
        # Write comment lines for headers
        csv_report.write("\n".join(["#" + s for s in headers])+"\n")
         # Write headers to the file
        csv_report.write("#dataset:" + "\n")
        # Write the output metric names
        csv_report.write("\t".join(s["output"]) + "\n")
    
    # Write data to the CSV report file                     
    csv_report.write("\t".join([str(s[field]) for field in s["output"]])+"\n")  # replace replace_nan by str
    csv_report.flush()      # Flush the buffer to ensure writing immediately
    nfiledone += 1          # Increment the count of processed files
    print('The results are updated.')
    

#############################################################################################################################
##########################                Print message box at the end of the code                 ##########################
#############################################################################################################################


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
##########################              Main function calling all the other functions              ##########################
#############################################################################################################################


if __name__ == '__main__':
    # Record the start time for runtime measurement
    start_time = time.time() 
    # Add the start time information to the headers list
    headers.append(f"start_time:\t{datetime.datetime.now()}")
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='')    # Path for the output folder name
    parser.add_argument('output_folder_name',
                        help = "the subfolder name on the '...\\UserInterface\\Data\\output_folder_name' directory.",
                        type=str)
    # Path for the input directory
    parser.add_argument('inputdir',
                        help = "input foldername consists of *.mha (*.nii or *.dcm) files. For example: 'E:\\Data\\Rectal\\input_data_folder'",
                        nargs = "*")
    parser.add_argument('-r', help="folders as name", default=False)
    parser.add_argument('-m', help="save foreground masks", default=False)
    parser.add_argument('-b', help="number of samples", default=1, type = int)
    parser.add_argument('-u', help="percent of middle images", default=100)
    parser.add_argument('-t', help="the address of the user-specified tags list (*.txt)", default=0)
    parser.add_argument('-c', help="if yes the ch computes objects", default=False)
    parser.add_argument('-s', help="type of scan (MRI or CT)", default='MRI', choices=['MRI', 'CT'])
    
    args = parser.parse_args() 
    root = args.inputdir[0]     # root = input_directory
    
    # Setting flags based on parsed arguments
    if args.r == 0:
        folders_flag = "False"
    else: 
        folders_flag = args.r
        
    if args.m == 0:
        save_masks_flag = "False" 
    else: 
        save_masks_flag = args.m
        
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
        
    if args.s == 0:
        scan_type = "MRI"
    else:
        scan_type = args.s
    
    # Les résultats vont dans le dossier UserInterface (dans MRQy)
    print_folder_note = os.getcwd() + os.sep + 'UserInterface'
    # Un dossier "Data" dans UserInterface va contenir les résultats de l'analyse des images
    fname_outdir = print_folder_note + os.sep + 'Data' + os.sep + args.output_folder_name
    
    overwrite_flag = "w"        
    headers.append(f"outdir:\t{os.path.realpath(fname_outdir)}")
    headers.append(f"scantype:\t{scan_type}")
    
    patients, names, dicom_spil, nondicom_spli, nondicom_names, mat_spli, mat_names = file_name(root)

    # Determine data types and set corresponding flags
    if len(dicom_spil) > 0 and len(nondicom_spli) > 0 and len(mat_spli) > 0:
        dicom_flag = True
        nondicom_flag = True
        mat_flag = True
    if len(dicom_spil) > 0 and len(nondicom_spli) > 0 and len(mat_spli) == 0:
        dicom_flag = True
        nondicom_flag = True
        mat_flag = False
    if len(dicom_spil) > 0 and len(nondicom_spli) == 0 and len(mat_spli) == 0:
        dicom_flag = True
        nondicom_flag = False
        mat_flag = False
    if len(dicom_spil) == 0 and len(nondicom_spli) > 0 and len(mat_spli) == 0:
        dicom_flag = False
        nondicom_flag = True
        mat_flag = False
    if len(dicom_spil) == 0 and len(nondicom_spli) > 0 and len(mat_spli) > 0:
        dicom_flag = False
        nondicom_flag = True
        mat_flag = True
    if len(dicom_spil) == 0 and len(nondicom_spli) == 0 and len(mat_spli) > 0:
        dicom_flag = False
        nondicom_flag = False
        mat_flag = True
    if len(dicom_spil) > 0 and len(nondicom_spli) == 0 and len(mat_spli) > 0:
        dicom_flag = True
        nondicom_flag = False
        mat_flag = True
    if len(dicom_spil) == 0 and len(nondicom_spli) == 0 and len(mat_spli) == 0:
        print('The input folder is empty or includes unsupported files format!')
    
    output_masks = os.path.join(fname_outdir, 'masks')
    os.makedirs(output_masks, exist_ok=True)

    
    # Process each data type (DICOM, non-DICOM, MAT) if available
    for i in range(len(names)):
        if dicom_flag:
            for j in range(len(dicom_spil)):
                v = volume_dicom(dicom_spil[j], names[j])
                folder_foregrounds = saveThumbnails_dicom(v, fname_outdir)
                s = BaseVolume_dicom(fname_outdir, v, j+1, folder_foregrounds, sample_size, ch_flag)
                worker_callback(s, fname_outdir)
            dicom_flag = False
            
        if nondicom_flag:
            for l,k in enumerate(nondicom_spli):
                mask = brain_extraction(k, output_masks)
                v = volume_nifti(k, nondicom_names[l])
                volume_masks = volume_nifti_masks(mask)
                saveThumbnails_nondicom(v, fname_outdir, volume_masks[0])
                s = BaseVolume_nondicom(fname_outdir, v, l+1, k, sample_size, ch_flag)
                worker_callback(s,fname_outdir)
            nondicom_flag = False
        
        if mat_flag:
            for j in range(len(mat_spli)):
                v = volume_mat(mat_spli[j], mat_names[j])
                folder_foregrounds = saveThumbnails_mat(v, fname_outdir)
                s = BaseVolume_mat(fname_outdir, v, j+1, folder_foregrounds)
                worker_callback(s, fname_outdir)
            mat_flag = False
    
    # Create the path for the result TSV file
    address = fname_outdir + os.sep + "results" + ".tsv" 
    
    # This creates the dataframe used for the tables, graphs, charts, UMAP and TSNE.
    # If you change things here (columns or rows used), it will modify what is shown in the user interface
    # Read the CSV file into a pandas dataframe, skipping the first two rows and using the third row as headers
    df = pd.read_csv(address, sep='\t', skiprows=4, header=0)
    df = df.drop(['Name of Images'], axis = 1)
    df = df.fillna('N/A')
    
    # Save processed data as IQM.csv
    df.to_csv(fname_outdir + os.sep + 'IQM.csv', index=False)
    print("The IQMs data are saved in the {} file. ".format(fname_outdir + os.sep + "IQM.csv"))
    
    # Print execution information
    print("Done!")
    print("MRQy program took", format((time.time() - start_time)/60, '.2f'), \
          "minutes for {} images to run.".format(len(names)))
    
    # Provide guidance for viewing the final results in the MRQy interface
    msg = "Please go to the '{}' directory and open up the 'index.html' file.\n".format(print_folder_note) + \
    "Click on 'View Results' and select '{}' file.\n".format(fname_outdir + os.sep + "results.tsv")   
    print_msg_box(msg, indent=3, width=None, title="To view the final MRQy interface results:")


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################