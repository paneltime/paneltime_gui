#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
sys.path.append(__file__.replace("__init__.py",''))
#import system_main as main

import gui
import gui_tempstore

def start():
	"""Starts the GUI"""
	gui_tempstore.test_and_repair()
	window=gui.window()
	window.mainloop() 
	
	
	
#Todo:
#(including paneltime)

#check that if works for no id and date variable
#add argument for null model (default: Y~Intercept)
#Put output functionality into the main_tabs object
#improve abort functionality, including what happens when tab is closed
#add durbin watson test:done
#output seems to be called twice
#change name of dataset (right click options?)
#make it possibel to add data by running 
#fix location issues with "+"-button. 
#create a right tab with a list of previoius estimates
#create right tab with all previously used and closed tabs available
#if one AR term is removed by reducing the AR order, the corespondig MA should be set to zero (if exists)
#have a backup for saved regressions and exe
#fix confusion about two option sets: one  in the starting environment and one in the data set
#Have a save symbol on all main tabs, so that the user can select a temporary folder
#save all files inside one zip-file
#check if Y is among the X variables
#add immediate command functionality to sub-pane
#add keyboard run shortcut and run selection
#make the dataset remember previoius alterations
#Add warning for un-nice hessian (avoid variables with huge variations in denomination)