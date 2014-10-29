# Cython header file of myfpzip.pyx

cdef extern from "fpzip.h":
    int fpzip_memory_write(void* buffer,int size,const void* data,const int * prec,int dp,int nx,int ny,int nz,int nf) 
    int fpzip_memory_read(const void* buffer,void* data,int * prec,int dp,int nx,int ny,int nz,int nf) 
