def Partition(varsg):
   return varsg


def Sync():
   return


def IsMaster():
   return True


def GetRank():
   return 0


def Sum(data):
   total = 0.0
   for k, v in data.iteritems():
      total = total + v
   return total


def Max(data):
   return data
