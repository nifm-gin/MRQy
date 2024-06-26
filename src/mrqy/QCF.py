import os
import numpy as np
import mrqy.QC as QC        # import QC
from scipy.signal import convolve2d as conv2
from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image
from skimage import exposure as ex
from skimage.filters import median
from skimage.morphology import square
from medpy.io import load
import warnings
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy

warnings.filterwarnings("ignore")


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

    def __init__(self, fname_outdir, v, ol, sample_size, ch_flag):
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
        # self.addToPrintList("ORIENTATION", orientation(scan), v, ol, 7)
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
        # self.addToPrintList("ORIENTATION", vol(v, sample_size, "ORIENTATION",fname_outdir, ch_flag), v, ol, 22)
        
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
        self.addToPrintList("CPP", vol(v, sample_size, "CPP", folder_foregrounds), v, ol, 8)
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
        F, B, c, f, b = foreground(I,outi_folder,v,i)
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
       

def foreground(img,save_folder,v,inumber):
    try:
        # Perform adaptive histogram equalization on the image
        h = ex.equalize_hist(img[:,:])*255
        # Binary thresholding using Otsu's method on the original and histogram-equalized images
        oi = np.zeros_like(img, dtype=np.uint16)            # Otsu thresholding of the original image
        oi[(img > threshold_otsu(img)) == True] = 1
        oh = np.zeros_like(img, dtype=np.uint16)            # Otsu trhesholding of the equalized image
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


################################################################
##############  Computing the different metrics  ###############
################################################################

# All the computed measurements are average values over the entire volume, which are calculated for every single slice separately.

# Mean of the foreground intensity values
def mean(F, B, c, f, b):
    return np.nanmean(f)

# Range of the foreground intensity values
def rang(F, B, c, f, b):
    return np.ptp(f)

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
    min_dim_index = image_shape.index(min(image_shape))
    return min_dim_index