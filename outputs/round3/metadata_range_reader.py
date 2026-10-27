"""Fill missing physical HDF5 metadata blocks on demand. Never read X logically."""
import io,hashlib,json
class MetadataRangeReader(io.RawIOBase):
 def __init__(self,path,url,fs,receipt,receipt_path):
  self.f=open(path,'r+b');self.url=url;self.fs=fs;self.receipt=receipt;self.rp=receipt_path;self.size=receipt['logical_file_size'];self.ranges=list(receipt['copied_segments'])+receipt.get('additional_metadata_blocks',[]);self.position=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.position
 def seek(self,offset,whence=0):
  self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset;return self.position
 def read(self,n=-1):
  n=self.size-self.position if n<0 else min(n,self.size-self.position);lo=self.position;hi=lo+n
  for begin in range((lo//1048576)*1048576,hi,1048576):
   end=min(begin+1048576,self.size)
   if any(begin>=a and end<=b for a,b in self.ranges):continue
   # Partial overlap still fetches one complete physical block, including
   # metadata chunks allocated before the obs object header.
   data=self.fs.cat_file(self.url,start=begin,end=end);assert len(data)==end-begin;self.f.seek(begin);self.f.write(data);self.f.flush();self.ranges.append([begin,end]);self.receipt.setdefault('additional_metadata_blocks',[]).append([begin,end]);self.receipt.setdefault('additional_metadata_block_hashes',[]).append({'start':begin,'end_exclusive':end,'sha256':hashlib.sha256(data).hexdigest()});self.rp.write_text(json.dumps(self.receipt,indent=2));print('METADATA_BLOCK_FILLED',begin,end,flush=True)
  self.f.seek(lo);data=self.f.read(n);self.position+=len(data);return data
 def readinto(self,b):
  data=self.read(len(b));b[:len(data)]=data;return len(data)
 def close(self):
  if not self.closed:self.f.close()
  super().close()
