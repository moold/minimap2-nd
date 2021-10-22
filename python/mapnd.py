import os
from ctypes import *

__all__ = ['Map']

class mm_idx_t (Structure):
	pass

class mm_out (Structure):
	_fields_ = [
		("hdr", c_char_p),
		("aln", c_char_p)
	]

MM = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "mapnd.so")
MM.build_index.argtypes = [c_char_p, c_char_p]
MM.build_index.restype = POINTER(mm_idx_t)
MM.destroy_index.argtypes = [POINTER(mm_idx_t)]
MM.map.argtypes = [POINTER(mm_idx_t), c_char_p, c_char_p]
MM.map.restype = POINTER(mm_out)
MM.destroy_out.argtypes = [POINTER(mm_out)]
MM.get_mid_occ.argtypes = [POINTER(mm_idx_t), c_char_p]
MM.get_mid_occ.restype = c_int32

class Map(object):
	def __init__(self, ref, options):
		self.options = options
		self.index = MM.build_index(ref.encode('utf-8'), options.encode('utf-8'))
		self.out = None

	def __del__(self):
		MM.destroy_index(self.index)

	def get_mid_occ(self):
		return MM.get_mid_occ(self.index, self.options.encode('utf-8'))

	def map(self, query, options=None):
		if not options:
			options = self.options
		out = MM.map(self.index, query.encode('utf-8'), options.encode('utf-8'))
		hdr = string_at(out.contents.hdr).decode("UTF-8") if out.contents.hdr else None
		aln = string_at(out.contents.aln).decode("UTF-8") if out.contents.aln else ''
		if hdr:
			self.out = '%s\n%s' % (hdr, aln)
		else:
			self.out = aln
		MM.destroy_out(out)
		return self.out.strip()

if __name__ == '__main__':
	options = "--sam-hit-only -Y -c -x map-hifi -t 10"
	ref = "human_g1k_v37.chr1.fasta"
	read = "chr1.hifi.fa"
	m = Map(ref, options)
	print (m.get_mid_occ())
	print (m.map(read, options))
