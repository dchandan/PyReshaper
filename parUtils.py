from mpi4py import MPI
def Partition(varsg):
   import math
   import numpy as np

   rank = MPI.COMM_WORLD.Get_rank()

   varsnames = []
   for k, v in varsg.iteritems():
      varsnames.append(k)

   glen = len(varsnames)
   nrank = MPI.COMM_WORLD.Get_size()

   llen = int(math.floor(glen / nrank))
   rem = glen - llen * nrank

   istart = np.arange(nrank)
   iend = np.arange(nrank)
   istart[0] = 0
   #print 'rem: ',rem
   #print 'llen: ',llen
   if 0 <= rem - 1:
      iend[0] = llen + 1
   else:
      iend[0] = llen

   for i in range(1, nrank):
      istart[i] = iend[i - 1] + 1 - 1
      if i <= rem - 1:
         iend[i] = istart[i] + llen + 1
      else:
         iend[i] = istart[i] + llen
   #print 'istart: ',istart
   #print 'iend: ',iend

   # Construct a local list of variables
   lvarnames = [ varsnames[i] for i in range(istart[rank], iend[rank])]

   # construct the local dictionary
   varsl = {}
   for k in lvarnames:
      varsl[k] = varsg[k]
   #print 'IAM: ',rank, 'variables: ',lvarnames

   return varsl


def Sync():
   MPI.COMM_WORLD.Barrier()


def IsMaster():
   rank = MPI.COMM_WORLD.Get_rank()
   if rank == 0:
      return True
   else:
      return False


def GetRank():
   return MPI.COMM_WORLD.Get_rank()


def Sum(data):
   total = 0.0
   for k, v in data.iteritems():
      total = total + v
   return MPI.COMM_WORLD.allreduce(total, op=MPI.SUM)


def Max(data):
   return MPI.COMM_WORLD.allreduce(data, op=MPI.MAX)
