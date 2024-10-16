##################### MRQy - Quality control and evaluation tool for MRI data #####################
###################################################################################################


###################################################################################################
######################################## Import libraries #########################################
###################################################################################################

import os
import numpy as np
import argparse
import datetime
import time
from medpy.io import load               # for .mha, .nii, or .nii.gz files
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pydicom                          # for .dcm files
from itertools import accumulate
import pandas as pd
import scipy
from scipy.cluster.vq import whiten
from scipy.signal import convolve2d as conv2
from sklearn.manifold import TSNE
from skimage import exposure as ex
from skimage.filters import threshold_otsu
from skimage.filters import median
from skimage.measure import find_contours
from skimage.morphology import convex_hull_image
from skimage.morphology import square
import math
import umap
from scipy.io import loadmat
import warnings        
warnings.filterwarnings("ignore")       # remove all warnings like conversion thumbnails

###################################################################################################
###################################################################################################

# Initialisation

nfiledone = 0
csv_report = None
first = True
headers = []


###################################################################################################
########################################### Definitions ###########################################
###################################################################################################


def patient_name(root):
    
    # Starting message
    print('MRQy is starting....')
    
    # Gathering relevant file paths based on extensions. MRQy supports .dcm, .mha, .nii, .gz and .mat files
    files = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(root) 
                for filename in filenames 
                if filename.endswith('.dcm') 
                or filename.endswith('.mha')
                or filename.endswith('.nii')
                or filename.endswith('.gz')
                or filename.endswith('.mat')]
    
    # Separating files based on their extensions
    mats = [i for i in files if i.endswith('.mat')]
    dicoms = [i for i in files if i.endswith('.dcm')]
    mhas = [i for i in files 
            if i.endswith('.mha')
            or i.endswith('.nii')
            or i.endswith('.gz')]
    
    # Extracting subject identifiers from the different files
    mhas_subjects = [os.path.basename(scan)[:os.path.basename(scan).index('.')] for scan in mhas]
    dicom_subjects = []
    mat_subjects = [os.path.basename(scan)[:os.path.basename(scan).index('.')] for scan in mats]
    
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
    subjects = subjects_id + mhas_subjects + mat_subjects
    # Displaying the total number of identified subjects
    print('The number of images is {}'.format(len(subjects)))
    # Returning various lists containing file paths, subjects, and DICOM splits
    return files, subjects, splits, mhas, mhas_subjects, mats, mat_subjects



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
             'ID': name_value,                                  # Patient / Subject ID
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



def volume_notdicom(scan, name):
    # Loading image data and header
    image_data, image_header = load(scan)
    # Get image dimensions
    image_shape = np.shape(image_data)
    # Get smallest dimension index : (a:b:c) with a --> sagittal, b --> coronal, c --> axial
    min_dim_index = image_shape.index(min(image_shape))
    # Extract the 2D images from the 3D image data
    '''
    if min_dim_index == 0:
        images = [image_data[i, :, :] for i in range(np.shape(image_data)[0])]      # Sagittal
    elif min_dim_index == 1:
        images = [image_data[:, i, :] for i in range(np.shape(image_data)[1])]      # Coronal
    elif min_dim_index == 2:
        images = [image_data[:, :, i] for i in range(np.shape(image_data)[2])]      # Axial
    else:
        images = [image_data[:, :, i] for i in range(np.shape(image_data)[2])]      # returns axial slices if it's a 3D volume
    '''
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



def volume_mat(mat_scan, name):
    # Loading volume data 
    v1 = loadmat(mat_scan)['vol']
    # Creating a dictionnary for the ID
    tags = {'ID': name}
    # Return the loaded volume data and the tags
    return v1, tags



def saveThumbnails_dicom(v, output):
    # Check if saving masks is enabled
    if save_masks_flag!='False':
        ffolder = output + '_foreground_masks'
        # Create a directory for foreground masks
        os.makedirs(ffolder + os.sep + v[1]['ID'])
    elif save_masks_flag=='False':
        ffolder = output
    # Create a directory for images
    os.makedirs(output + os.sep + v[1]['ID'])       # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    # Save images as thumbnails
    for i in range(0, len(v[0]), sample_size):
        plt.imsave(output + os.sep + v[1]['ID'] + os.sep + v[1]['ID'] + '(%d).png' % i, v[0][i], cmap = cm.Greys_r)                             # Le chemin + nom qu'ont les images.png dans leur dossier respectif
    # Print the number of saved images and the directory
    print('The number of %d images are saved to %s' % (len(v[0]),output + os.sep + v[1]['ID']))
    return ffolder + os.sep + v[1]['ID']



def saveThumbnails_mat(v, output):
    # Check if saving masks is enabled
    if save_masks_flag!='False':
        ffolder = output + '_foreground_masks'
        # Create a directory for foreground masks
        os.makedirs(ffolder + os.sep + v[1]['ID'])
    elif save_masks_flag=='False':
        ffolder = output
    # Create a directory for images 
    os.makedirs(output + os.sep + v[1]['ID'])       # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    # Save image as thumbnails
    for i in range(np.shape(v[0])[2]):
        plt.imsave(output + os.sep + v[1]['ID']+ os.sep + v[1]['ID'] + '(%d).png' % int(i+1), v[0][:,:,i], cmap = cm.Greys_r)                   # Le chemin + nom qu'ont les images.png dans leur dossier respectif
        # Print the number of saved images and the directory
    print('The number of %d images are saved to %s' % (np.shape(v[0])[2],output + os.sep + v[1]['ID']))
    return ffolder + os.sep + v[1]['ID']



def saveThumbnails_nondicom(v, output):
    # Create a directory for images
    os.makedirs(output + os.sep + v[1])             # Dans le dossier "Data", il va y avoir la création des dossiers contenant les images au format png
    # Save images as thumbnails, with a rotation
    for i in range(len(v[0])):
        # Save the images. Because of the precedent image processing operations, the image was rotated 90° in anti-clockwise. So we apply a 90° clockwise rotation
        plt.imsave(output + os.sep + v[1] + os.sep + v[1] + '(%d).png' % int(i+1), scipy.ndimage.rotate(v[0][i],90), cmap = cm.Greys_r)         # Le chemin + nom qu'ont les images.png dans leur dossier respectif
    # print('image number %d out of %d is saved to %s' % (int(i+1), len(v[0]),output + os.sep + v[1]))
    print('The number of %d images are saved to %s' % (len(v[0]),output + os.sep + v[1]))



class BaseVolume_dicom(dict):

    def __init__(self, fname_outdir, v, ol, folder_foregrounds, sample_size, ch_flag):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        # Add patient information to the output list
        self.addToPrintList("Patient", v[1]['ID'], v, ol, 170)
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
        # Display information about the patient's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the patient with the name of <%s> is %s' % (ol,il,name, v[1]['ID'], val))



class BaseVolume_nondicom(dict):

    def __init__(self, fname_outdir, v, ol, scan, sample_size, ch_flag):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        
        # Add information to the output list
        self.addToPrintList("Patient", v[1], v, ol, 170)
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
        self.addToPrintList("MEAN", vol(v, sample_size, "Mean", fname_outdir, ch_flag), v, ol, 8)
        self.addToPrintList("RNG", vol(v, sample_size, "Range", fname_outdir, ch_flag), v, ol, 9)
        self.addToPrintList("VAR", vol(v, sample_size, "Variance", fname_outdir, ch_flag), v, ol, 10)
        self.addToPrintList("CV", vol(v, sample_size, "CV", fname_outdir, ch_flag), v, ol, 11)
        self.addToPrintList("CPP", vol(v, sample_size, "CPP", fname_outdir, ch_flag), v, ol, 12)
        self.addToPrintList("PSNR", vol(v, sample_size, "PSNR", fname_outdir, ch_flag), v, ol, 13)
        self.addToPrintList("SNR1", vol(v, sample_size, "SNR1", fname_outdir, ch_flag), v, ol, 14)
        self.addToPrintList("SNR2", vol(v, sample_size, "SNR2", fname_outdir, ch_flag), v, ol, 15)
        self.addToPrintList("SNR3", vol(v, sample_size, "SNR3", fname_outdir, ch_flag), v, ol, 16)
        self.addToPrintList("SNR4", vol(v, sample_size, "SNR4", fname_outdir, ch_flag), v, ol, 17)
        self.addToPrintList("CNR", vol(v, sample_size, "CNR", fname_outdir, ch_flag), v, ol, 18)
        self.addToPrintList("CVP", vol(v, sample_size, "CVP", fname_outdir, ch_flag), v, ol, 19)
        self.addToPrintList("CJV", vol(v, sample_size, "CJV", fname_outdir, ch_flag), v, ol, 20)
        self.addToPrintList("EFC", vol(v, sample_size, "EFC", fname_outdir, ch_flag), v, ol, 21)
        self.addToPrintList("FBER", vol(v, sample_size, "FBER", fname_outdir, ch_flag), v, ol, 22)
        
    def addToPrintList(self, name, val, v, ol, il):
        # Add a new key-value pair to the dictionary
        self[name] = val
        self["output"].append(name)
        # Display information about the patient's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the patient with the name of <%s> is %s' % (ol,il,name, v[1], val))



class BaseVolume_mat(dict):

    def __init__(self, fname_outdir, v, ol,folder_foregrounds, sample_size):
        # Initialize the dictionary
        dict.__init__(self)

        # Initialize attributes within the dictionary
        self["warnings"] = [] 
        self["output"] = []
        
        # Add patient information to the output list
        self.addToPrintList("Patient", v[1]['ID'], v, ol, 170)
        
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
        # Display information about the patient's metrics
        if name != 'Name of Images' and il != 170:
            print('%s-%s. The %s of the patient with the name of <%s> is %s' % (ol,il,name, v[1]['ID'], val))



def vol(v, sample_size, kk, outi_folder, ch_flag):
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
    for i in range(1, len(v[0]), sample_size):
        I = v[0][i]
        # I = I - np.min(I)  # for CT 
        # Calculate foreground and background intensities
        F, B, c, f, b = foreground(I, outi_folder, v, i)
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



def foreground(img, save_folder, v, inumber):
    try:
        # Perform adaptive histogram equalization on the image
        h = ex.equalize_hist(img[:,:])*255
        # Binary thresholding using Otsu's method on the original and histogram-equalized images
        oi = np.zeros_like(img, dtype=np.uint16)            # Otsu thresholding of the original image
        oi[(img > threshold_otsu(img)) == True] = 1
        oh = np.zeros_like(img, dtype=np.uint16)            # Otsu thresholding of the equalized image
        oh[(h > threshold_otsu(h)) == True] = 1
        # Compute weights for the images based on the thresholding results
        nm = img.shape[0] * img.shape[1]
        w1 = np.sum(oi)/(nm)                # w1 = weight of Otsu on original image
        w2 = np.sum(oh)/(nm)                # w2 = weight of Otsu on equalized image
        # Compute a new image using the calculated weights
        ots = np.zeros_like(img, dtype=np.uint16)
        new = (w1 * img) + (w2 * h)
        ots[(new > threshold_otsu(new)) == True] = 1 
        # Obtain convex hull of the thresholded image
        conv_hull = convex_hull_image(ots)
        # Calculate the foreground and background images based on the convex hull
        ch = np.multiply(conv_hull, 1)
        fore_image = ch * img
        back_image = (1 - ch) * img
    except Exception:
        # If an exception occurs, return default values
        fore_image = img.copy()
        back_image = np.zeros_like(img, dtype=np.uint16)
        conv_hull = np.zeros_like(img, dtype=np.uint16)
        ch = np.multiply(conv_hull, 1)
        
    # if not os.path.isdir(save_folder + os.sep + v[1]['ID']):                              
    return fore_image, back_image, conv_hull, img[conv_hull], img[conv_hull==False]         # fore_image = F    back_image = B      conv_hull = c       img[conv_hull] = f      img[conv_hull==False] = b



# All the computed measurements are average values over the entire volume, which are calculated for every single slice separately.

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
        csv_report.write("#dataset:"+"\t".join(s["output"])+"\n")
    
    # Write data to the CSV report file                     
    csv_report.write("\t".join([str(s[field]) for field in s["output"]])+"\n")          # replace replace_nan by str
    csv_report.flush()      # Flush the buffer to ensure writing immediately
    nfiledone += 1          # Increment the count of processed files
    print('The results are updated.')



def tsv_to_dataframe(tsvfileaddress):       # Creates a dataframe with all the data extracted from metadata and calculated from images
    # This creates the dataframe used for the tables, graphs, charts, UMAP and TSNE.
    # If you change things here (columns or rows used), it will modify what is shown in the user interface
    return pd.read_csv(tsvfileaddress, sep='\t', skiprows=2, header=0)          # Read the CSV file into a pandas dataframe, skipping the first two rows and using the third row as headers



def data_whitening(dframe):
    # Fill missing values with N/A in the dataframe
    dframe = dframe.fillna('N/A')
    # Create a copy of the dataframe     
    df = dframe.copy()
    # Select columns excluding object type
    df = df.select_dtypes(exclude=['object'])
    # Enter the names of the columns you want included for the UMAP and TSNE 
    df = df[["MEAN", "RNG", "VAR", "CV", "CPP", "PSNR","SNR1", "SNR2", "SNR3", "SNR4", "CNR", "CVP", "CJV", "EFC", "FBER"]]
    # Apply whitening transformation to the selected numeric data. The data calculated from the images are transformed to be used by UMAP and TSNE
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



###################################################################################################
########################## Main function calling all the other functions ##########################
###################################################################################################


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
    parser.add_argument('-s', help="save foreground masks", default=False)
    parser.add_argument('-b', help="number of samples", default=1, type = int)
    parser.add_argument('-u', help="percent of middle images", default=100)
    parser.add_argument('-t', help="the address of the user-specified tags list (*.txt)", default=0)
    parser.add_argument('-c', help="if yes the ch computes objects", default=False)
    
    args = parser.parse_args() 
    root = args.inputdir[0]
    
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
    
    # print(os.getcwd())
    print_folder_note = os.getcwd() + os.sep + 'UserInterface'                                  # Les résultats vont dans le dossier UserInterface (dans MRQy)
    # print(print_folder_note)
    fname_outdir = print_folder_note + os.sep + 'Data' + os.sep + args.output_folder_name       # Un dossier "Data" dans UserInterface va contenir les résultats de l'analyse des images
    
    overwrite_flag = "w"        
    headers.append(f"outdir:\t{os.path.realpath(fname_outdir)}") 
    patients, names, dicom_spil, nondicom_spli, nondicom_names, mat_spli, mat_names = patient_name(root)

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
    
    # Process each data type (DICOM, non-DICOM, MAT) if available
    for i in range(len(names)):
        if dicom_flag:
            for j in range(len(dicom_spil)):
                v = volume_dicom(dicom_spil[j], names[j])
                folder_foregrounds = saveThumbnails_dicom(v,fname_outdir)
                s = BaseVolume_dicom(fname_outdir, v,j+1,folder_foregrounds, sample_size, ch_flag)
                worker_callback(s,fname_outdir)
            dicom_flag = False
            
        if nondicom_flag:
            for l,k in enumerate(nondicom_spli):
                v = volume_notdicom(k, nondicom_names[l])
                saveThumbnails_nondicom(v,fname_outdir)
                s = BaseVolume_nondicom(fname_outdir, v, l+1, k, sample_size, ch_flag)
                worker_callback(s,fname_outdir)
            nondicom_flag = False
        
        if mat_flag:
            for j in range(len(mat_spli)):
                v = volume_mat(mat_spli[j], mat_names[j])
                folder_foregrounds = saveThumbnails_mat(v,fname_outdir)
                s = BaseVolume_mat(fname_outdir, v,j+1,folder_foregrounds)
                worker_callback(s,fname_outdir)
            mat_flag = False
    
    # Create the path for the result TSV file
    address = fname_outdir + os.sep + "results" + ".tsv" 
    

    if len(names) < 6:
        # t-SNE and UMPA cannot be performed if we have less than 6 images
        print('Skipped the t-SNE and UMAP computation because of insufficient data. The UMAP and t-SNE process need at least 6 input data.')
        df = tsv_to_dataframe(address)
    else:    
        df = cleanup(address, 30)
        df = df.drop(['Name of Images'], axis=1)
        df = df.rename(columns={"#dataset:Patient": "Patient", "x":"TSNEX","y":"TSNEY", "u":"UMAPX", "v":"UMAPY" })
        df = df.fillna('N/A')

    # Save processed data as IQM.csv
    df.to_csv(fname_outdir + os.sep + 'IQM.csv', index=False)
    print("The IQMs data are saved in the {} file. ".format(fname_outdir + os.sep + "IQM.csv"))
    
    # Draw the contours of the ROI after the analysis is done so that the contours don't interfere with the measures
    # This is useful if you want to verify that that the ROI includes the head/brain/bodypart so that the metrics are calculated correctly
    def image_contour(input_folder, output_folder):
        print("Drawing the contours of the ROI... Almost done")
        # Verify the files in the directory
        for root, dirs, files in os.walk(input_folder):
            output_root = os.path.join(output_folder, os.path.relpath(root, input_folder))
            os.makedirs(output_root, exist_ok=True)
            for file in files:
                if file.endswith('csv') or file.endswith('.tsv'):
                    continue
                input_path = os.path.join(root, file)
                output_path = os.path.join(output_root, file)
                img = plt.imread(input_path)
                # Get input image shape
                input_height, input_width = img.shape[:2]
                # Convert to grayscale
                if img.shape[2] == 4:
                    img = img[:,:,:3]   # exclude the alpha channel
                img = img.mean(axis=2)
                # Perform adaptive equalization of the image 
                ae = ex.equalize_adapthist(img)
                # oi = Otsu thresholding of the original image
                oi = np.zeros_like(img, dtype=np.uint16)            
                oi[(img > threshold_otsu(img)) == True] = 1
                # oa = Otsu thresholding of the adaptive equalized image
                oa = np.zeros_like(img, dtype=np.uint16)            
                oa[(ae > threshold_otsu(ae)) == True] = 1
                # Compute weights for the images nased on the thresholding results
                nm = img.shape[0] * img.shape[1]
                w1 = np.sum(oi)/(nm)        # w1 = weight of Otsu on original image (img)
                w2 = np.sum(oa)/(nm)        # w2 = weight of Otsu on adaptive equalized (ae)
                # New image using the weights of the Otsu thresholding and adaptive equalization
                ots = np.zeros_like(img, dtype=np.uint16)                     
                new = (w1 * img) + (w2 * ae)
                ots[(new > threshold_otsu(new)) == True] = 1
                # Convex hull of the new image using the weights of the Otsu thresholding and adaptive equalization (ots_2)
                conv_hull = convex_hull_image(ots)
                # Compute foreground image
                ch = np.multiply(conv_hull, 1)
                fore_image = ch * img
                # Make sure the output image has the same shape as input image
                fore_image = fore_image[:input_height, :input_width]
                # Defintion to draw the contours
                def contours(conv_hull):
                    # Find contours of the foreground image (= the convex_hull)
                    contours = find_contours(conv_hull)
                    for contour in contours:
                        plt.plot(contour[:, 1], contour[:, 0], linewidth=3, color='green')
                        plt.axis('off')
                # contours(fore_image, conv_hull)
                plt.imshow(fore_image, cmap='gray') ; contours(conv_hull)
                plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
                plt.close()
    
    input_folder = fname_outdir                             # Where the .png are saved
    output_folder = fname_outdir                            # Where the _contours.png are saved
    image_contour(input_folder, output_folder)
    
    
    # Print execution information
    print("Done!")
    print("MRQy program took", format((time.time() - start_time)/60, '.2f'), \
          "minutes for {} subjects and the overal {} MRI slices to run.".format(len(names),len(patients)))
    
    # Provide guidance for viewing the final results in the MRQy interface
    msg = "Please go to the '{}' directory and open up the 'index.html' file.\n".format(print_folder_note) + \
    "Click on 'View Results' and select '{}' file.\n".format(fname_outdir + os.sep + "results.tsv")   
    print_msg_box(msg, indent=3, width=None, title="To view the final MRQy interface results:")