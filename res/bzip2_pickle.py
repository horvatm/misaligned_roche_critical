import pickle
import bz2

def save(object, filename, protocol = pickle.HIGHEST_PROTOCOL):
        """Saves a compressed object to disk
        """
        file = bz2.BZ2File(filename, 'wb')
        file.write(pickle.dumps(object, protocol))
        file.close()

def load(filename):
        """Loads a compressed object from disk
        """
        file = bz2.BZ2File(filename, 'rb')
        buffer = ""
        while True:
                data = file.read()
                if data == "":
                        break
                buffer += data
        object = pickle.loads(buffer)
        file.close()
        return object


def load2dict(filename):
  """
    Loading data into a dictionary
  """
  #
  # Read all data
  #
  data = np.loadtxt(filename)

  qs = np.unique(data[:,0]) # set of all considered q
  Fs = np.unique(data[:,1])  # set of all considered F
  ths = np.unique(data[:,2]) #  # set of all considered theta

  #~ print "q=", qs
  #~ print "F=", Fs
  #~ print "th=", ths

  #
  # Slice all data into chucks of constant (q,F)
  #

  ch = defaultdict(list)
  for d in data:
    q, F = d[:2]
    ch[(q,F)].append(d[-4:]) 

  #
  # Compactify chucks and save
  #
  for k, v in ch.items():
    ch[k] = np.array(v)
  
  return ch

def save_dict(ch, filename):
  compress.save(ch, filename)


def load_dict(filename):
  return compress.load(filename)

