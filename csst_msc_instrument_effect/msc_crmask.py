#!/usr/bin/env python3

##############################################################
#Author  Hu Yi (NAOC, huyi.naoc@gmail.com), Zhang Yajie (TJU)#
##############################################################

#Thanks to Zhang Keming who is the author of deepCR.
#Thanks to the authors of astropy and ccdproc.
#Software package dependencies, numpy, astropy, ccdproc and deepCR
#Installation dependencies on Ubuntu 20.04
#apt install python3-numpy python3-astropy python3-ccdproc python3-pip
#python3 -m pip install deepCR
#Version 0.1

import configparser
import numpy as np
from astropy.io import fits as pyfits
from astropy import units as u
from astropy.nddata import CCDData
from pathlib import Path
import ccdproc 

from deepCR import deepCR
from deepCR import train
import torch

import datetime
import sys


__all__ = ['CRMask']

class CRMask:
    def __init__(self, obj, flag = None, mask = None, sky = None, save_flag = True, update_flag = True, save_name = None, flag_suffix = 'flg', clean_suffix = 'crclean', append_flag = False, mask_suffix = 'crmask', model = 'deepCR', fill_flag = True, fill_method = 'inpainting', gpu_flag = False, config_path = 'MSC_crmask.ini', **kwargs):
        """
            Instantiation of CRMask with specified model configuration.
            Parameters
            ----------
            obj : string, Path, astropy.io.fits.HDUList, numpy.ndarray, astropy.nddata.CCDData or list of string
                if model is ``deepCR``, ``lacosmic``, obj is input image to be cosmic-ray masked
                if model is ``deepCR_train``, obj is input training and validating images
            flag : (optional) string, Path, astropy.io.fits.HDUList, numpy.ndarray or astropy.nddata.CCDData
                flag image(s), default is None
            mask : (optional) string or list of string
                mask image(s), default is None, necessary when model is ``deepCR_train``, otherwise, let it alone
            sky : (optional) string or list of string
                sky image(s), default is None, optional when model is ``deepCR_train``, otherwise, let it alone 
            save_flag : (optional) boolean
                whether save CR mask (and cleaned) image, default is True
            update_flag : (optional) boolean
                whether update flag image, default is True
            save_name : (optional) string
                output mask, cleaned and flag filename. default is None. If save_name is None, use the filename as the input. And if save_flag is True, save_name is None, and obj is a numpy.ndarray or a astropy.nddata.CCDData, cr_mask will raise a ValueError exception
            flag_suffix : (optional) string
                suffix name of flag file, if flag is a numpy.ndarray or a astropy.nddata.CCDData, default is ``flg``
            mask_suffix : (optional) string
                suffix name of mask file, default is ``crmask``
            clean_suffix : (optional) string
                suffix name of cleaned file, default is ``crclean``
            model : (optional) string
                model type, can be ``deepCR``, ``lacosmic``, ``deepCR_train``, default is ``deepCR``
            fill_flag : (optional) boolean
                whehter generate cleaned image, default is True
            fill_method : (optional) string
                fill method for CR contaminated pixel, can be ``inpainting``, ``meanmask``, ``meanmed``, default is ``inpainting``
            gpu_flag : (optional) boolean
                whether use GPU, default is False
            config_path : (optional) string
                configuration file path, default is ``../conf/MSC_crmask.ini``
        """
        if model == 'deepCR_train':
            self.image_sets = obj
            self.mask_sets = mask
            self.ignore_sets = flag
            self.sky_sets = sky
            self.gpu_flag = gpu_flag
            if config_path != None:
                config = configparser.ConfigParser()
                config.read(config_path)
                if config.has_option('global', 'gpu_flag'):
                    self.append_flag = config.getboolean('global', 'gpu_flag')
                self.config = config
            return

        if isinstance(obj, str) or isinstance(obj, Path):
            self.hdulist = pyfits.open(obj)
        elif isinstance(obj, pyfits.HDUList):
            self.hdulist = obj
        elif isinstance(obj, np.ndarray) or isinstance(obj, CCDData):
            hdulist = pyfits.HDUList()
            primary = pyfits.PrimaryHDU()
            header = primary.header
            header['NEXTEND'] = 1
            hdulist.append(primary)
            image = pyfits.ImageHDU()
            header = image.header
            header['BITPIX'] = -32
            header['NAXIS'] = 2
            header.insert('NAXIS', ('NAXIS1', obj.shape[0]), after = True)
            header.insert('NAXIS1', ('NAXIS2', obj.shape[1]), after = True)
            header['EXTNAME'] = 'SCI'
            header['EXTVER'] = 1
            if isinstance(obj, CCDData):
                image.data = np.asarray(CCDData)
            else:
                image.data = obj
            hdulist.append(image)
            self.hdulist = hdulist
        else:
            raise TypeError('For cosmic ray masking and cleaning, obj must be a string, or a pathlib.Path object, or a astropy.io.fits.HDUList object!')  
        if flag != None:
            if isinstance(flag, str) or isinstance(flag, Path):
                flag_hdulist = pyfits.open(flag, mode = 'update')
            elif isinstance(flag, pyfits.HDUList):
                flag_hdulist = flag
            elif isinstance(flag, np.ndarray) or isinstance(flag, CCDData):
                flag_hdulist = pyfits.HDUList()
                flag_primary = pyfits.PrimaryHDU(header = self.hdulist[0].header)
                flag_hdulist.append(flag_hdulist)
                if isinstance(obj, CCDData):
                    flag_data = np.asarray(CCDData)
                    flag_image = pyfits.ImageHDU(data=flag_data, header=self.hdulist[1].header)
                else:
                    flag_image = pyfits.ImageHDU(data=flag, header=self.hdulist[1].header)
                flag_hdulist.append(flag_image)
            else:
                raise TypeError('For cosmic ray masking and cleaning, mask must be a string, or a pathlib.Path object, or a astropy.io.fits.HDUList object!')
            self.flag = flag_hdulist
        else:
            self.flag = None

        self.update_flag = update_flag
        self.save_flag = save_flag
        self.save_name = save_name
        self.flag_suffix = flag_suffix
        self.clean_suffix = clean_suffix
        self.append_flag = append_flag
        self.mask_suffix = mask_suffix
        self.gpu_flag = gpu_flag
        self.fill_method = fill_method
        self.model = model

        if config_path != None:
            config = configparser.ConfigParser()
            config.read(config_path)

            if config.has_option('global', 'save_flag'):
                self.save_flag = config.getboolean('global', 'save_flag')

            if config.has_option('global', 'clean_suffix'):
                self.clean_suffix = config.get('global', 'clean_suffix')

            if config.has_option('global', 'append_flag'):
                self.append_flag = config.getboolean('global', 'append_flag')

            if config.has_option('global', 'flag_suffix'):
                self.flag_suffix = config.get('global', 'flag_suffix')
            
            if config.has_option('global', 'mask_suffix'):
                self.mask_suffix = config.get('global', 'mask_suffix')

            if config.has_option('global', 'gpu_flag'):
                self.append_flag = config.getboolean('global', 'gpu_flag')
            
            if config.has_option('global', 'fill_flag'):
                self.fill_flag = config.getboolean('global', 'fill_flag')
            
            if config.has_option('global', 'fill_method'):
                self.fill_method = config.get('global', 'fill_method')

            if config.has_option('global', 'model'):
                self.model = config.get('global', 'model')

            if config.has_option('global', 'torch_thread'):
                self.torch_thread = config.get('global', 'torch_thread')
                torch.set_num_threads(self.torch_thread)
            self.config = config

    def cr_mask_lacosmic(self):
        
        config = self.config
        
        if config.has_option('lacosmic', 'sigclip'):
            sigclip = config.getfloat('lacosmic', 'sigclip')
        else:
            sigclip = 4.0

        if config.has_option('lacosmic', 'sigfrac'):
            sigfrac = config.getfloat('lacosmic', 'sigfrac')
        else:
            sigfrac = 0.3
        
        if config.has_option('lacosmic', 'objlim'):
            objlim = config.getfloat('lacosmic', 'objlim')
        else:
            objlim = 5.0
        
        if config.has_option('lacosmic', 'gain'):
            gain = config.getfloat('lacosmic', 'gain')
        else:
            gain = 1.0

        if config.has_option('lacosmic', 'readnoise'):
            readnoise = config.getfloat('lacosmic', 'readnoise')
        else:
            readnoise = 6.5

        if config.has_option('lacosmic', 'satlevel'):
            satlevel = config.getfloat('lacosmic', 'satlevel')
        else:
            satlevel = 65535.0

        if config.has_option('lacosmic', 'pssl'):
            pssl = config.getfloat('lacosmic', 'pssl')
        else:
            pssl = 0.0
    
        if config.has_option('lacosmic', 'niter'):
            niter = config.getint('lacosmic', 'niter')
        else:
            niter = 4
        
        if config.has_option('lacosmic', 'sepmed'):
            sepmed = config.getboolean('lacosmic', 'sepmed')
        else:
            sepmed = True


        if config.has_option('lacosmic', 'cleantype'):
            cleantype = config.get('lacosmic', 'cleantype')
        else:
            cleantype = self.fill_method

        if config.has_option('lacosmic', 'fsmode'):
            fsmode = config.get('lacosmic', 'fsmode')
        else:
            fsmode = 'median'

        if config.has_option('lacosmic', 'psfmodel'):
            psfmodel = config.get('lacosmic', 'psfmodel')
        else:
            psfmodel = 'gauss'

        if config.has_option('lacosmic', 'psffwhm'):
            psffwhm = config.getfloat('lacosmic', 'psffwhm')
        else:
            psffwhm = 2.5

        if config.has_option('lacosmic', 'psfsize'):
            psfsize = config.getint('lacosmic', 'psfsize')
        else:
            psfsize = 7

        if config.has_option('lacosmic', 'psfk'):
            psfk = config.get('lacosmic', 'psfk')
        else:
            psfk = None

        if config.has_option('lacosmic', 'psfbeta'):
            psfbeta = config.getfloat('lacosmic', 'psfbeta')
        else:
            psfbeta = 4.765

        if config.has_option('lacosmic', 'gain_apply'):
            gain_apply = config.getboolean('lacosmic', 'gain_apply')
        else:
            gain_apply = True

        data = self.hdulist[1].data

        start_time = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
        cleaned, masked = ccdproc.cosmicray_lacosmic(data, sigclip, sigfrac, objlim, gain, readnoise, satlevel, pssl, niter, sepmed, cleantype, fsmode, psfmodel, psffwhm, psfsize, psfk, psfbeta, False, gain_apply)
        end_time = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
        masked_hdulist = pyfits.HDUList()
        masked_primary = pyfits.PrimaryHDU(header=self.hdulist[0].header)
        masked_hdulist.append(masked_primary)
        masked = masked.astype(np.uint8)
        masked_image = pyfits.ImageHDU(data=masked, header=self.hdulist[1].header)
        masked_hdulist.append(masked_image)
        
        #TODO
        #add history keywords here
        mask = None
        if self.flag != None:
            self.flag[1].data |= (masked<<4)
            self.flag[1].header.add_history('Use ccdproc.cosimicray_lacosmicfor cosmic ray detecting')
            value = 'CRMask start at {0}'.format(start_time)
            self.flag[1].header.add_history(value)
            value = 'CRMask end   at {0}'.format(end_time)
            self.flag[1].header.add_history(value)
        
        masked_hdulist[1].header.add_history('Use ccdproc.cosimicray_lacosmic for cosmic ray detecting')
        value = 'CRMask start at {0}'.format(start_time)
        masked_hdulist[1].header.add_history(value)
        value = 'CRMask end at   {0}'.format(end_time)
        masked_hdulist[1].header.add_history(value)

        if self.fill_flag:
            cleaned_hdulist = pyfits.HDUList()
            cleaned_primary = pyfits.PrimaryHDU(header=self.hdulist[0].header)
            cleaned_hdulist.append(cleaned_primary)
            cleaned_image = pyfits.ImageHDU(data=cleaned, header=self.hdulist[1].header)
            cleaned_hdulist.append(cleaned_image)
            cleaned_hdulist[1].header.add_history('Use ccdproc.cosimicray_lacosmic for cosmic ray detecting')
            value = 'CRMask start at {0}'.format(start_time)
            cleaned_hdulist[1].header.add_history(value)
            value = 'CRMask end at   {0}'.format(end_time)
            cleaned_hdulist[1].header.add_history(value)

            return masked_hdulist, cleaned_hdulist
        else:
            return masked_hdulist

    def cr_mask_deepCR(self):

        config = self.config
        
        if config.has_option('deepCR', 'threshold'):
            threshold = config.getfloat('deepCR', 'threshold')
        else:
            threshold = 0.5

        if config.has_option('deepCR', 'inpaint'):
            inpaint = config.getboolean('deepCR', 'inpaint')
        else:
            inpaint = True

        if config.has_option('deepCR', 'binary'):
            binary = config.getboolean('deepCR', 'binary')
        else:
            binary = True

        if config.has_option('deepCR', 'patch'):
            patch = config.getint('deepCR', 'patch')
        else:
            patch = 256

        if config.has_option('deepCR', 'segment'):
            segment = config.getboolean('deepCR', 'segment')
        else:
            segment = True
        
        if config.has_option('deepCR', 'parallel'):
            parallel = config.getboolean('deepCR', 'parallel')
        else:
            parallel = True

        if config.has_option('deepCR', 'clean_model'):
            clean_model = config.get('deepCR', 'clean_model')
        else:
            clean_model = 'ACS-WFC-F606W-2-32'
        
        if config.has_option('deepCR', 'inpaint_model'):
            inpaint_model = config.get('deepCR', 'inpaint_model')
        else:
            inpaint_model = 'ACS-WFC-F606W-2-32'
        
        if config.has_option('deepCR', 'n_jobs'):
            n_jobs = config.getint('deepCR', 'n_jobs')
        else:
            n_jobs = -1

        if self.gpu_flag:
            model = deepCR(clean_model, inpaint_model, device = 'GPU', hidden=50)
        else:
            model = deepCR(clean_model, inpaint_model, device = 'CPU', hidden=50)
        
        data = self.hdulist[1].data

        masked_hdulist = pyfits.HDUList()
        masked_primary = pyfits.PrimaryHDU(header=self.hdulist[0].header)
        masked_hdulist.append(masked_primary)
        if inpaint:
            cleaned_hdulist = pyfits.HDUList()
            cleaned_primary = pyfits.PrimaryHDU(header=self.hdulist[0].header)
            cleaned_hdulist.append(cleaned_primary)

        start_time = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
        if binary:
            if inpaint:
                masked, cleaned = model.clean(data, threshold = threshold, inpaint = inpaint, segment = segment, patch = patch, parallel = parallel, n_jobs = n_jobs)
            else:
                masked = model.clean(data, threshold = threshold, inpaint = inpaint, segment = segment, patch = patch, parallel = parallel, n_jobs = n_jobs)
        else:
            if inpaint:
                masked, cleaned = model.clean(data, threshold = threshold, inpaint = inpaint, binary = False, segment = segment, patch = patch, parallel = parallel, n_jobs = n_jobs)
            else:
                masked = model.clean(data, threshold = threshold, inpaint = inpaint, binary = False, segment = segment, patch = patch, parallel = parallel, n_jobs = n_jobs)
        end_time = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f") 
        
        masked =  masked.astype(np.uint8)
        masked_image = pyfits.ImageHDU(data=masked, header=self.hdulist[1].header)
        masked_hdulist.append(masked_image)
        
        #TODO
        #add history keywords here
        if self.flag != None:
            self.flag[1].data |= (masked<<4)
            self.flag[1].header.add_history('Use deepCR for cosmic ray detecting')
            value = 'CRMask start at {0}'.format(start_time)
            self.flag[1].header.add_history(value)
            value = 'CRMask end at   {0}'.format(end_time)
            self.flag[1].header.add_history(value)
        masked_hdulist[1].header.add_history('Use deepCR for cosmic ray detecting')
        value = 'CRMask start at {0}'.format(start_time)
        masked_hdulist[1].header.add_history(value)
        value = 'CRMask end at   {0}'.format(end_time)
        masked_hdulist[1].header.add_history(value)
        if inpaint:
            cleaned_image = pyfits.ImageHDU(data=cleaned, header=self.hdulist[1].header)
            cleaned_hdulist.append(cleaned_image)
            cleaned_hdulist[1].header.add_history('Use deepCR for cosmic ray detecting')
            value = 'CRMask start at {0}'.format(start_time)
            cleaned_hdulist[1].header.add_history(value)
            value = 'CRMask end   at {0}'.format(end_time)
            cleaned_hdulist[1].header.add_history(value)
            return masked_hdulist,cleaned_hdulist
        else:
            return masked_hdulist

    #return a hdulist of a masked image, and a hdulist of a cleaned image if fill_flag is set.
    def cr_mask(self):
        #here do cr mask task
        if self.model  == 'lacosmic':
            if self.fill_flag:
                masked, cleaned = CRMask.cr_mask_lacosmic(self)
            else:
                masked = CRMask.cr_mask_lacosmic(self)
        elif self.model == 'deepCR':
            if self.fill_flag:
                masked, cleaned = CRMask.cr_mask_deepCR(self)
            else:
                masked = CRMask.cr_mask_deepCR(self)

        else:
            raise ValueError('Cosmic ray model are not supported!')

        #save the result to a file.
        output = None
        if self.save_flag:
            if self.save_name != None:
                output = self.save_name
            else:
                output = self.hdulist.filename()
            #append to the input
            if self.append_flag:
                self.hdulist.append(masked[1])
                if self.fill_flag:
                    self.hdulist.append(cleaned[1])
                if output != None:
                    self.hdulist.writeto(output)
                else:
                    #TODO
                    #raise an exception.
                    pass
            else:
                if output != None:
                    if output.split('.')[-1] == 'fits':
                        mask_output = output[0:-4] + self.mask_suffix + '.fits'
                        clean_output = output[0:-4] + self.clean_suffix + '.fits'
                    else:
                        mask_output = output + '.' + self.mask_suffix + '.fits'
                        clean_output = output + '.' + self.clean_suffix + '.fits'
                    masked.writeto(mask_output)
                    if self.fill_flag:
                        cleaned.writeto(clean_output)
                else:
                    #TODO
                    #raise an exception here
                   pass

        #update flag file.
        if self.update_flag and self.flag != None:
            flag_output = self.flag.filename()
            if flag_output == None:
                if output != None:
                    if output.split('.')[-1] == 'fits':
                        flag_output = output[0:-4] + self.flag_suffix + '.fits'
                    else:
                        flag_output = output + '.' + self.flag_suffix + '.fits'
                else:
                    #raise an exception here.
                    pass
            self.flag.writeto(flag_output)

        if self.fill_flag:
            return (masked, cleaned)
        else:
            return masked

    def cr_train_deepCR_image_to_ndarray(self, image_sets, patch):
        
        if isinstance(image_sets, str):
            if image_sets[:-4] == '.npy':
                input_image = np.load(image_sets)
            elif image_sets[:-4] == '.fits':
                image = pyfits.open(image_sets)
                data = image[1].data
                x, y = data.shape
                m = x // patch
                n = y // patch
                start_x = np.mod(x, patch) // 2
                start_y = np.mod(y, patch) // 2
                input_image = data[start_x:start_x + m * patch][ start_y:start_y + m * patch].reshape[m*n][patch][patch]
            else:
                #TODO
                #raise an exception
                pass
        if isinstance(image_sets, list):    
            input_list = []
            for image_file in image_sets:
                if isinstance(image_sets, str):
                    if image_file[:-4] == '.npy':
                        input_image = np.load(image_file)
                    elif image_file[:-4] == '.fits':
                        image = pyfits.open(image_file)
                        data = image[1].data
                        x, y = data.shape
                        m = x // patch
                        n = y // patch
                        start_x = np.mod(x, patch) // 2
                        start_y = np.mod(y, patch) // 2
                        input_image = data[start_x:start_x + m * patch][start_y:start_y + m * patch].reshape[m*n][patch][patch]
                    else:
                        #TODO
                        #raise an exception
                        pass
                    input_list.append(input_image)
                else:
                    #TODO
                    #raise an exception
                    pass
            input_image = np.concatenate(input_list)
        return input_image

    def cr_train_deepCR_prepare_data(self, patch):
        if self.image_sets != None:
            self.training_image = CRMask.cr_train_deepCR_image_to_ndarray(self, self.image_sets, patch)
        else:
            raise ValueError('training image should not be None')

        if self.mask_sets != None:
            self.training_mask = CRMask.cr_train_deepCR_image_to_ndarray(self, self.mask_sets, patch)
        else:
            raise ValueError('mask image should not be None')
        
        if self.flag != None:
            self.training_ignore = CRMask.cr_train_deepCR_image_to_ndarray(self, self.flag, patch)
        else:
            self.training_ignore = None

        if self.sky != None:
            self.training_sky = CRMask.cr_train_deepCR_image_to_ndarray(self, self.sky, patch)
        else:
            self.training_sky = None

    
    def cr_train_deepCR(self):
        config = self.config

        if config.has_option('deepCR', 'aug_sky'):
            aug_sky = np.array(config.get('deepCR', 'aug_sky').split()).astype(np.float32)
            aug_sky0 = aug_sky[0]
            aug_sky1 = aug_sky[1]
        else:
            aug_sky0 = 0.
            aug_sky1 = 0.

        if  config.has_option('deepCR', 'model_name'):
            model_name = config.get('deepCR', 'model_name')
        else:
            model_name = 'model'

        if  config.has_option('deepCR', 'hidden'):
            hidden = config.getint('deepCR', 'hidden')
        else:
            hidden = 50
            
        if  config.has_option('deepCR', 'epoch'):
            epoch = config.getint('deepCR', 'epoch')
        else:
            epoch = 50

        if config.has_option('deepCR', 'patch'):
            patch = config.getint('deepCR', 'patch')
        else:
            patch = 256

        if config.has_option('deepCR', 'batch'):
            batch = config.getint('deepCR', 'batch')
        else:
            batch = 16

        if config.has_option('deepCR', 'lr'):
            lr = config.getfloat('deepCR', 'lr')
        else:
            lr = 0.005

        if config.has_option('deepCR', 'auto_lr_decay'):
            auto_lr_decay = config.getboolean('deepCR', 'auto_lr_decay')
        else:
            auto_lr_decay = True

        if config.has_option('deepCR', 'lr_decay_patience'):
            lr_decay_patience = config.getint('deepCR', 'lr_decay_patience')
        else:
            lr_decay_patience = 4
        
        if config.has_option('deepCR', 'lr_decay_factor'):
            lr_decay_factor = config.getfloat('deepCR', 'lr_decay_factor')
        else:
            lr_decay_factor = 0.1

        if config.has_option('deepCR', 'save_after'):
            save_after = config.getfloat('deepCR', 'save_after')
        else:
            save_after = 100000.0

        if config.has_option('deepCR', 'plot_every'):
            plot_every = config.getint('deepCR', 'plot_every')
        else:
            plot_every = 4

        if config.has_option('deepCR', 'verbose'):
            verbose = config.getboolean('deepCR', 'verbose')
        else:
            verbose = True
        
        if config.has_option('deepCR', 'use_tqdm'):
            use_tqdm = config.getboolean('deepCR', 'use_tqdm')
        else:
            use_tqdm = False

        if config.has_option('deepCR', 'use_tqdm_notebook'):
            use_tqdm_notebook = config.getboolean('deepCR', 'use_tqdm_notebook')
        else:
            use_tqdm_notebook = False
        
        if  config.has_option('deepCR', 'directory'):
            directory = config.get('deepCR', 'directory')
        else:
            directory = './'
        
        CRMask.cr_train_deepCR_prepare_data(self, patch)
        trainer = train(self.training_image, self.training_mask, self.training_ignore, self.training_sky, [aug_sky0, aug_sky1], model_name, hidden, self.gpu_flag, epoch, batch, lr, auto_lr_decay, lr_decay_patience, lr_decay_patience, save_after, plot_every, verbose, use_tqdm, use_tqdm_notebook, directory)
        trainer.train()

    def cr_train(self):
        if self.model == 'deepCR_train':
            CRMask.cr_train_deepCR(self)
        else:
            raise ValueError('Unsupported training model')

#to run the code
#./MSC_crmask xxxx.fits deepCR for deepCR
#./MSC_crmask xxxx.fits lacosmic for lacosmic

if __name__ == '__main__':
    crobj = CRMask(sys.argv[1], model = sys.argv[2])
    crobj.cr_mask()
