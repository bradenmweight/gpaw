import os
import sys
import time
from hdf5_highlevel import File, HyperslabSelection

import numpy as np

intsize = 4
floatsize = np.array([1], float).itemsize
complexsize = np.array([1], complex).itemsize
itemsizes = {'int': intsize, 'float': floatsize, 'complex': complexsize}

from gpaw.mpi import broadcast

class Writer:
    def __init__(self, name, comm=None):
        self.comm = comm # not used much in writer
        self.dims = {}        
        try:
           if self.comm.rank == 0:
               if os.path.isfile(name):
                   os.rename(name, name[:-5] + '.old'+name[-5:])
           self.comm.barrier()
        except AttributeError:
           if os.path.isfile(name):
               os.rename(name, name[:-5] + '.old'+name[-5:])

        self.file = File(name, 'w', self.comm.get_c_object())
        self.dims_grp = self.file.create_group("Dimensions")
        self.params_grp = self.file.create_group("Parameters")
        self.file.attrs['title'] = 'gpaw_io version="0.1"'
        
    def dimension(self, name, value):
        if name in self.dims.keys() and self.dims[name] != value:
            raise Warning('Dimension %s changed from %s to %s' % \
                          (name, self.dims[name], value))
        self.dims[name] = value
        self.dims_grp.attrs[name] = value

    def __setitem__(self, name, value):
        # if self.comm.rank == 0: writing on master leads to hang
        self.params_grp.attrs[name] = value

    def add(self, name, shape, array=None, dtype=None, 
            parallel=False, write=True):
        if array is not None:
            array = np.asarray(array)

        # self.dtype, type, itemsize = self.get_data_type(array, dtype)
        if dtype is None:
            self.dtype = array.dtype
        else:
            self.dtype = dtype

        shape = [self.dims[dim] for dim in shape]
        if not shape:
            shape = [1,]
        self.dset = self.file.create_dataset(name, shape, self.dtype)
        if array is not None:
            self.fill(array, parallel=parallel, write=write)

    def fill(self, array, *indices, **kwargs):

        parallel = kwargs.pop('parallel', False)
        write = kwargs.pop('write', True)

        if parallel:
            collective = True
        else:
            collective = False

        if not write:
            selection = None
        elif indices: 
            selection = HyperslabSelection(indices, self.dset.shape)
        else:
            selection = 'all'
        self.dset.write(array, selection, collective)            

    def get_data_type(self, array=None, dtype=None):
        if dtype is None:
            dtype = array.dtype

        if dtype in [int, bool]:
            dtype = np.int32

        dtype = np.dtype(dtype)
        type = {np.int32: 'int',
                np.float64: 'float',
                np.complex128: 'complex'}[dtype.type]

        return dtype, type, dtype.itemsize

    def append(self, name):
        raise NotImplementedError('Append with HDF5 not available.')

    def close(self):
        mtime = int(time.time())
        self.file.attrs['mtime'] = mtime
        self.dims_grp.close()
        self.params_grp.close()
        self.dset.close()
        self.file.close()
        
class Reader:
    def __init__(self, name, comm=None):
        self.comm = comm # used for broadcasting replicated data 
        self.file = File(name, 'r', self.comm.get_c_object())
        self.params_grp = self.file['Parameters']
        self.hdf5_reader = True

    def dimension(self, name):
        dims_grp = self.file['Dimensions']
        return dims_grp.attrs[name]
    
    def __getitem__(self, name):
        obj = None
        if self.comm.rank == 0: # not sure that this is necessary
            obj = self.params_grp.attrs[name]
              
        value = broadcast(obj, 0, self.comm)

        try:
            value = eval(value, {})
        except (SyntaxError, NameError, TypeError):
            pass
        return value

    def has_array(self, name):
        return name in self.file.keys()
    
    def get(self, name, *indices, **kwargs):

        parallel = kwargs.pop('parallel', False)
        read = kwargs.pop('read', True)
        out = kwargs.pop('out', None)
        assert not kwargs

        if parallel:
            collective = True
        else: 
            collective = False

        dset = self.file[name]
        if indices:
            selection = HyperslabSelection(indices, dset.shape)
            mshape = selection.mshape
        else:
            selection = 'all'
            mshape = dset.shape

        if not read:
            selection = None

        if out is None:
            array = np.ndarray(mshape, dset.dtype, order='C')
        else:
            assert type(out) is np.ndarray
            # XXX Check the shapes are compatible
            assert out.shape == mshape
            assert out.dtype == dset.dtype
            array = out

        dset.read(array, selection, collective)

        if array.shape == ():
            return array.item()
        else:
            return array

    def get_reference(self, name, *indices):
        dset = self.file[name]
        array = dset[indices]
        return array

    def get_parameters(self):
        return self.params_grp.attrs
    
    def close(self):
        self.file.close()
