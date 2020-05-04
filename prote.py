from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from RegexParser import *
import pandas as pd
from bs4 import BeautifulSoup
import numpy as np
import json
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
import re
import random
import time




def main(url,name):
 

	driver = webdriver.Chrome(ChromeDriverManager().install())
	# Open the URL
	driver.get(url)
	time.sleep(1)  # Let the user actually see something!
	content = driver.page_source

	# getting the html content of this url
	soup = BeautifulSoup(content, 'html.parser')

	# gene = parameter.split('/')[3].split('-')[1]
	# gene_dic={}
	# celline_dic={}
	# block_matches = my_parser.parse("<g class=\"bar_g\"(.*?)<br>Category",driver.page_source)
	# count=0
	# l=[]

	# find input box one
	print(gene_name[i])
	searchbox_One = driver.find_element_by_id("react-select-2-input")
	searchbox_One.send_keys('Gene')
	searchbox_One.send_keys(Keys.ENTER)
	print('selected gene')

	time.sleep(1)
	searchbox_Two = driver.find_element_by_id("react-select-3-input")
	searchbox_Two.send_keys(f'{name}')
	time.sleep(4)
	searchbox_Two.send_keys(Keys.RETURN)
	print('entered gene')

	time.sleep(2)
	searchbox_Three = driver.find_element_by_id("react-select-4-input")
	searchbox_Three.send_keys('CRISPR (Avana) Public 20Q1')
	searchbox_Three.send_keys(Keys.ENTER)
	print('entered dataset')

	time.sleep(1)
	# show precomputed associations
	checkbox = driver.find_element_by_name("associationTable")
	checkbox.click()
	time.sleep(2)

	# unclick dataset

	print('unclicking dataset')

	checkboxOne = driver.find_element_by_name("Copy Number (Absolute)")
	checkboxOne.click()
	checkboxTwo = driver.find_element_by_name("Copy Number Public 20Q1")
	checkboxTwo.click()
	checkboxThree = driver.find_element_by_name("Expression Public 20Q1")
	checkboxThree.click()
	checkboFour = driver.find_element_by_name("Mutation Public 20Q1")
	checkboFour.click()
	# Get data
	time.sleep(1)
	divs = driver.find_elements_by_class_name('bottom-row-oriented')

	for div in divs:
		button = div.find_element_by_css_selector('button.btn.btn-default.btn-sm')
		button.click()
		print('downloaded')

	time.sleep(2)
	driver.close()


if __name__ == "__main__":

	gene_name=[]
	# list of gene name to want to fetch
	with open("/Users/anna/Desktop/Animalcrossing/gene_list.txt", "r") as f:
		for line in f:
		 gene_name.append(line.strip())

	print(gene_name)


		
	url = 'https://depmap.org/portal/interactive/'

	my_parser = RegexParser(url)
	f = open("tt.txt", "w")

	for i in range(len(gene_name)):

		try:
			
			main(url, gene_name[i])
			
		except Exception:
			print("Oops!  That was no web")
		
