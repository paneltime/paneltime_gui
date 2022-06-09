#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy import stats as scstats
from paneltimegui import gui_functions as guif
import os
import shutil
from paneltime import stat_functions as stat
from matplotlib import pyplot  as plt




class process_charts(ttk.Frame):
	def __init__(self,window,master,main_tabs,tabs):
		style = ttk.Style()
		style.configure("TFrame", background='white')		
		ttk.Frame.__init__(self,master,style='new.TFrame')
		self.window=window
		self.subplot=tabs.subplot
		self.print_subplot=tabs.print_subplot
		self.generate = GenerateCharts()
		self.add_content()

		
	def add_content(self):
		self.n_charts=3
		self.columnconfigure(0,weight=1)
		for i in range(self.n_charts+1):
			self.rowconfigure(i,weight=1)
			
		tk.Label(self,text='Charts on normalized residuals:',bg='white',font='Tahoma 10 bold').grid(row=0,column=0)			

		self.charts=[]
		for i in range(self.n_charts):
			frm=tk.Frame(self,background='white')
			frm.rowconfigure(0,weight=1)
			frm.rowconfigure(1)
			frm.columnconfigure(0,weight=1)
			self.charts.append(tk.Label(frm,background='white'))
			self.charts[i].grid(row=0,column=0)	
			chart_path=os.path.join(os.getcwd(),'img',f'chart{i}.png')
			self.charts[i].path=chart_path
			guif.setbutton(frm, 'Save image', lambda: self.save(self.n_charts-i-1),bg='white').grid(row=1,column=0)
			frm.grid(row=i+1)
		
	def save(self,i,f=None):
		if not hasattr(self.charts[i],'path'):
			print('No graphics displayed yet')
			return
		name=self.charts[i].name
		if f is None:
			f = tk.filedialog.asksaveasfile(mode='bw', defaultextension=".jpg",initialfile=f"{name}.jpg")		
		if f is None:
			return
		shutil.copyfile(self.charts[i].path, f)
		
	def plot(self,panel, ll):
		self.generate.save_all(panel, ll)
		for c in self.charts:
			plot_to_chart(c)

def plot_to_chart(chart):
	return#currently disabled, in development
	if hasattr(chart,'graph_file'):
		chart.graph_file.close()
	chart.graph_file=Image.open(chart.path)
	img = ImageTk.PhotoImage(chart.graph_file,master=chart)
	chart.configure(image=img)
	chart.graph_img=img	




#!/usr/bin/env python
# -*- coding: utf-8 -*-





class GenerateCharts():
	def __init__(self):
		return#currently disabled, in development
		self.subplot=plt.subplots(1,figsize=(4,2.5),dpi=75)
		self.chart_list=[
			['histogram',self.histogram],
			['correlogram',self.correlogram],
			['correlogram_variance',self.correlogram_variance]
		]

	def save_all(self, panel, ll):
		return#currently disabled, in development
		for name,chart in self.chart_list:
			chart(panel, ll,self.subplot,f'img/{name}.png')				

	def histogram(self, panel,ll,subplot,f):
		return#currently disabled, in development
		N,T,k = panel.X.shape
		fgr,axs=subplot
		n = ll.e_norm_centered.shape[2]
		e = ll.e_norm_centered[panel.included[2]].flatten()
		N = e.shape[0]
		e = e.reshape((N,1))

		grid_range = 4
		grid_step = 0.05	
		h,grid = histogram(e,grid_range,grid_step)
		norm = scstats.norm.pdf(grid)*grid_step	

		axs.bar(grid,h,color='grey', width=0.025,label='histogram')
		axs.plot(grid,norm,'green',label='normal distribution')
		axs.legend(prop={'size': 6})
		name='Histogram - frequency'
		axs.set_title(name)
		save(subplot,f)

	def correlogram(self, panel,ll,subplot,f):
		return#currently disabled, in development
		fgr,axs=subplot
		lags=20
		rho=stat.correlogram(panel, ll.e_norm_centered,lags)
		x=np.arange(lags+1)
		axs.bar(x,rho,color='grey', width=0.5,label='correlogram')
		name='Correlogram - residuals'
		axs.set_title(name)
		save(subplot,f)

	def correlogram_variance(self, panel,ll,subplot,f):
		return#currently disabled, in development
		N,T,k=panel.X.shape
		fgr,axs=subplot
		lags=20
		e2=ll.e_norm_centered**2
		e2=(e2-anel.mean(e2))*panel.included[3]
		rho=stat.correlogram(panel, e2,lags)
		x=np.arange(lags+1)
		axs.bar(x,rho,color='grey', width=0.5,label='correlogram')
		name='Correlogram - squared residuals'
		axs.set_title(name)
		save(subplot,f)

def histogram(x,grid_range,grid_step):
	N,k=x.shape
	grid_n=int(2*grid_range/grid_step)
	grid=np.array([i*grid_step-grid_range for i in range(grid_n)]).reshape((1,grid_n))
	ones=np.ones((N,1))
	x_u=np.concatenate((ones,x>=grid),1)
	x_l=np.concatenate((x<grid,ones),1)
	grid=np.concatenate((grid.flatten(),[grid[0,-1]+grid_step]))
	histogram=np.sum((x_u*x_l),0)
	if int(np.sum(histogram))!=N:
		raise RuntimeError('Error in histogram calculation')
	return histogram/N,grid



def save(subplot,save_file):
	fgr,axs=subplot
	fgr.savefig(save_file)
	axs.clear()