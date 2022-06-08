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