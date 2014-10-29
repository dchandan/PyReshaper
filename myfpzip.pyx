cimport myfpzip
cimport numpy
cimport cython
import numpy


# Python wrapper for fpzip_memory_write and fpzip_memory_read

# Call fpzip compress/decompress for 1d float data
def my_python_fpzip_1d(numpy.ndarray[float,ndim=1,mode="c"] input_array not None,precision_bit,dp,input_shape):
    
    cdef numpy.ndarray[numpy.float32_t,ndim=1,mode="c"] x
    cdef numpy.ndarray[numpy.float32_t,ndim=1,mode="c"] y
    cdef numpy.ndarray[numpy.float32_t,ndim=1,mode="c"] z
    buffer_array=numpy.empty(input_shape,dtype=float,order='C')
    output_array=numpy.empty(input_shape,dtype=float,order='C')
    cdef int size=1,nf,pre[1]
    cdef numpy.ndarray ndim =numpy.ones((5,),dtype=int)
    j=0
    for i in input_shape:
       ndim[j]=i
       j=j+1
       size=size*i
    pre[0]=precision_bit
    x=numpy.ascontiguousarray(input_array,dtype=numpy.float32)
    y=numpy.ascontiguousarray(buffer_array,dtype=numpy.float32)
    z=numpy.ascontiguousarray(output_array,dtype=numpy.float32)
    
    result=myfpzip.fpzip_memory_write(&y[0],size*4,&x[0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    result=myfpzip.fpzip_memory_read(&y[0],&z[0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    return result,z

# Call fpzip compress/decompress for 2d float data
def my_python_fpzip_2d(numpy.ndarray[float,ndim=2,mode="c"] input_array not None,precision_bit,dp,input_shape):
    
    cdef numpy.ndarray[numpy.float32_t,ndim=2,mode="c"] x
    cdef numpy.ndarray[numpy.float32_t,ndim=2,mode="c"] y
    cdef numpy.ndarray[numpy.float32_t,ndim=2,mode="c"] z
    buffer_array=numpy.empty(input_shape,dtype=float,order='C')
    output_array=numpy.empty(input_shape,dtype=float,order='C')
    cdef int size=1,nf,pre[1]
    cdef numpy.ndarray ndim =numpy.ones((5,),dtype=int)
    j=0
    for i in input_shape:
       ndim[j]=i
       j=j+1
       size=size*i
    pre[0]=precision_bit
    x=numpy.ascontiguousarray(input_array,dtype=numpy.float32)
    y=numpy.ascontiguousarray(buffer_array,dtype=numpy.float32)
    z=numpy.ascontiguousarray(output_array,dtype=numpy.float32)
    
    result=myfpzip.fpzip_memory_write(&y[0,0],size*4,&x[0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    result=myfpzip.fpzip_memory_read(&y[0,0],&z[0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    return result,z

# Call fpzip compress/decompress for 3d float data
def my_python_fpzip_3d(numpy.ndarray[float,ndim=3,mode="c"] input_array not None,precision_bit,dp,input_shape):
    
    cdef numpy.ndarray[numpy.float32_t,ndim=3,mode="c"] x
    cdef numpy.ndarray[numpy.float32_t,ndim=3,mode="c"] y
    cdef numpy.ndarray[numpy.float32_t,ndim=3,mode="c"] z
    buffer_array=numpy.empty(input_shape,dtype=float,order='C')
    output_array=numpy.empty(input_shape,dtype=float,order='C')
    cdef int size=1,nf,pre[1]
    cdef numpy.ndarray ndim =numpy.ones((5,),dtype=int)
    j=0
    for i in input_shape:
       ndim[j]=i
       j=j+1
       size=size*i
    pre[0]=precision_bit
    x=numpy.ascontiguousarray(input_array,dtype=numpy.float32)
    y=numpy.ascontiguousarray(buffer_array,dtype=numpy.float32)
    z=numpy.ascontiguousarray(output_array,dtype=numpy.float32)
    
    result=myfpzip.fpzip_memory_write(&y[0,0,0],size*4,&x[0,0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    result=myfpzip.fpzip_memory_read(&y[0,0,0],&z[0,0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    return result,z


# Call fpzip compress/decompress for 3d double data
def my_python_fpzip_3d_double(numpy.ndarray[float,ndim=3,mode="c"] input_array not None,precision_bit,dp,input_shape):
    
    cdef numpy.ndarray[numpy.float64_t,ndim=3,mode="c"] x
    cdef numpy.ndarray[numpy.float64_t,ndim=3,mode="c"] y
    cdef numpy.ndarray[numpy.float64_t,ndim=3,mode="c"] z
    buffer_array=numpy.empty(input_shape,dtype=float,order='C')
    output_array=numpy.empty(input_shape,dtype=float,order='C')
    cdef int size=1,nf,pre[1]
    cdef numpy.ndarray ndim =numpy.ones((5,),dtype=int)
    j=0
    for i in input_shape:
       ndim[j]=i
       j=j+1
       size=size*i
    pre[0]=precision_bit
    x=input_array.view(dtype=numpy.float64)
    y=buffer_array.view(dtype=numpy.float64)
    z=output_array.view(dtype=numpy.float64)
    
    result=myfpzip.fpzip_memory_write(&y[0,0,0],size*4,&x[0,0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    result=myfpzip.fpzip_memory_read(&y[0,0,0],&z[0,0,0],pre,dp,ndim[0],ndim[1],ndim[2],1)
    return result,z

# can be called from fastjuice.py to start fpzip compression/decompression
def my_python_fpzip(input_array,precision_bit,dp):
    if input_array.ndim == 1:
      com_size,recon_data=my_python_fpzip_1d(input_array,precision_bit,dp,input_array.shape)
    elif input_array.ndim == 2:
      com_size,recon_data=my_python_fpzip_2d(input_array,precision_bit,dp,input_array.shape)
    elif input_array.ndim == 3:
      com_size,recon_data=my_python_fpzip_3d(input_array,precision_bit,dp,input_array.shape)
    return com_size,recon_data
