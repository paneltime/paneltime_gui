#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pickle
import tempfile
import os
import zipfile

tdr=tempfile.gettempdir()


fname_window=os.path.join(tdr,'paneltime.win')
fname_datasets=os.path.join(tdr,'paneltime.datasets')

max_sessions=20
file_name_list=[fname_window, fname_datasets]



def load_obj(fname):
	for i in [0,1]:
		try:
			f=open(fname, "r+b")
			u= pickle.Unpickler(f)
			u=u.load()
			f.close()
			return u 
		except Exception as e:
			print(e)
			recreate_from_zip()
			if i==1:
				return
	
def save_zip():
	wdr=os.getcwd()
	zip_file_path=os.path.join(wdr,'data.paneltime')
	zip_arch=zipfile.ZipFile(zip_file_path,'w')	
	for f in file_name_list:
		if os.path.isfile(f):
			zip_arch.write(f,os.path.basename(f))
	zip_arch.close()
	
def test_and_repair():
	ok=True
	for f in file_name_list:
		ok=ok and os.path.isfile(f)
		if not ok:
			break
	if not ok:
		recreate_from_zip()
	
def recreate_from_zip():
	wdr=os.getcwd()
	zip_file_path=os.path.join(wdr,'data.paneltime')
	try:
		zip_arch=zipfile.ZipFile(zip_file_path,'r')
	except:
		return
	for f in zip_arch.filelist:
		zf=zip_arch.read(f.filename)
		fl=open(os.path.join(tdr,f.filename),'wb')
		fl.write(zf)
		fl.close()
	zip_arch.close()
	
def save_obj(fname,obj):
	f=open(fname, "w+b")
	pickle.dump(obj,f)   
	f.flush() 
	f.close()	





				
			
		
	