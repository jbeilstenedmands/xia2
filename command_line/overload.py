from __future__ import division

def read_cbf_image(cbf_image):
  from cbflib_adaptbx import uncompress
  import binascii

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(cbf_image, 'rb').read()
  data_offset = data.find(start_tag) + 4
  cbf_header = data[:data_offset - 4]

  fast = 0
  slow = 0
  length = 0

  for record in cbf_header.split('\n'):
    if 'X-Binary-Size-Fastest-Dimension' in record:
      fast = int(record.split()[-1])
    elif 'X-Binary-Size-Second-Dimension' in record:
      slow = int(record.split()[-1])
    elif 'X-Binary-Number-of-Elements' in record:
      length = int(record.split()[-1])
    elif 'X-Binary-Size:' in record:
      size = int(record.split()[-1])

  assert(length == fast * slow)

  pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                            fast = fast, slow = slow)

  return pixel_values

def get_overload(cbf_file):
  for record in open(cbf_file):
    if 'Count_cutoff' in record:
      return int(record.split()[-2])

def build_hist():
  import sys
  from scitbx.array_family import flex

  hist = { }

  thresh = 100
  scale = thresh / get_overload(sys.argv[1])

  for image in sys.argv[1:]:
    sys.stdout.write('.')
    sys.stdout.flush()
    data = read_cbf_image(image)
    scaled = (scale * data.as_double()).iround()
    for d in scaled:
      if not d in hist:
        hist[d] = 0
      hist[d] += 1

  print '\n'
  for d in sorted(hist)[1:]:
    if d < 0.1 * thresh:
      print d, hist[d]
    else:
      print d, hist[d], '***'

if __name__ == '__main__':
  build_hist()
