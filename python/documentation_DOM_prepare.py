#!/usr/bin/env python

from __future__ import print_function
import sys, os
import copy
from bs4 import BeautifulSoup
from collections import OrderedDict

import regex as re

from argparse import ArgumentParser, RawDescriptionHelpFormatter

def set(d, set, value):
    keys = list(set)
    keys.append(value)
    #latest = keys.pop()
    for k in keys:
        d = d.setdefault(k, OrderedDict())
        
# <span id="document-Introduction"></span><div class="section" id="introduction">
# <h2>Introduction<a class="headerlink" href="#introduction" title="Permalink to this headline"></a></h2>
# <p>Slamdunk maps and analyzes SLAM-Seq data.</p>
# </div>
# <span id="document-Installation"></span><div class="section" id="installation">
# <h2>Installation<a class="headerlink" href="#installation" title="Permalink to this headline"></a></h2>
# <p>This section covers the installation of <em>slamdunk</em>. Alternatively one could also run <em>slamdunk</em> from out of the box <a class="reference internal" href="index.html#docker-label"><em>Docker</em></a> containers.</p>
# <p>There are 2 different possibilities:</p>
# <ol class="arabic simple">
# <li>Installation from <a class="reference external" href="https://pypi.python.org/pypi">PyPI</a> using the <a class="reference internal" href="index.html#pip-label"><em>Python Package Index</em></a> <strong>(recommended)</strong></li>
# <li>Installation from <a class="reference internal" href="index.html#source-label"><em>Source</em></a></li>
# </ol>
# <div class="section" id="requirements">
# <h3>Requirements<a class="headerlink" href="#requirements" title="Permalink to this headline"></a></h3>
# <p>There are no major requirements for <em>slamdunk</em>. The python package will acquire all external dependencies by itself.</p>
# <div class="section" id="r-runtime">

# index.html#document-Introduction
# index.html#document-Installation
#            index.html#requirements
#                index.html#r-runtime
#            index.html#python-package-index

# <div class="docs_section">
# <h1 class="section-header" id="document-Introduction"><a href="#document-Introduction" class="header-link"><span class="glyphicon glyphicon-link"></span></a>Introduction</h1>
# <div class="docs_block" id="Introduction.rst"><p>Slamdunk maps and analyzes SLAM-Seq data.</p>
# </div>
# </div>
        
        
def doStuff(soup, section, TOCMap):
    
    body = soup.find('div',class_="toctree-wrapper compound")
    
    tag = body.find(id=section)

    if tag.name == "span":
         
# <span id="document-Introduction"></span><div class="section" id="introduction">
# <h2>Introduction<a class="headerlink" href="#introduction" title="Permalink to this headline"></a></h2>
         
        content = tag.next_sibling
        content['class'] = "docs_section"
        del content['id']
        sectionName = content.h2.contents[0]
        header = content.h2
        header.name = "h1"
        header['class'] = "section-header"
        header['id'] = section
        header.contents[0] = ""
        permalink = content.a
        permalink['class'] = "header-link"
        permalink['href'] = "#" + section
        del permalink['title']
        permalink.string = ""
        
        glyphlink = soup.new_tag("span")
        glyphlink["class"] = "glyphicon glyphicon-link"
        permalink.append(glyphlink)
        
        header.append(str(sectionName))
        
        data = soup.new_tag("div", class_="docs_block", id = section)
        
        sibling = glyphlink.parent.parent.next_sibling
        
        while sibling :
            data.append(sibling)
            sibling = glyphlink.parent.parent.next_sibling
            
        content.append(data)
        
        tag.decompose()
        
        # FROM
        
#     <div class="section" id="requirements">
#     <h3>Requirements<a class="headerlink" href="#requirements" title="Permalink to this headline"></a></h3>
        # TO
#     <h1 id="requirements"><a href="#requirements" class="header-link"><span class="glyphicon glyphicon-link"></span></a>Requirements</h1>
#     <div class="section" id="requirements">

    elif tag.name == "div":
        
        header = tag.find(re.compile("^h"))
        header.name = "h" + str(TOCMap[section] - 1)
        header["id"] = section
        sectionName = header.contents[0]
        
        permalink = header.a
        permalink['class'] = "header-link"
        permalink['href'] = "#" + section
        del permalink['title']
        permalink.string = ""
        
        header.contents[0] = ""
        
        glyphlink = soup.new_tag("span")
        glyphlink["class"] = "glyphicon glyphicon-link"
        permalink.append(glyphlink)
        
        header.append(str(sectionName))
        
        print(header)
        sys.stdin.readline()
        
        
        
        
            #new_ol.insert(0, data.extract())

        #content.append(docs_block_wrapper)
        #print(docs_block_wrapper)
        #print(str(soup)[1:3000])
        #sys.stdin.readline()     
         
#         html += "<div class=\"docs_section\">"
#         html += "<h1 class=\"" + section + "\" id=\"document-Introduction\">"
#         html += "<a href=\"#" + section + "\" class=\"header-link\">"
#         html += "<a href=\"#" + section + "\" class=\"header-link\">"
#        html += "<span class=\"glyphicon glyphicon-link\"></span></a>" + header +"</h1>"
#        html += "<div class=\"docs_block\">"
    
def walkthroughTOC(soup, toc, TOCMap, level = 0):
    for k, v in toc.iteritems():
        doStuff(soup, k, TOCMap)
        #print('\t' * level, end="")
        #print(k)
        if v:
            walkthroughTOC(soup, v, TOCMap, level + 1)


 # Info
usage = "Parsing SingleHtml sphinx-build into documentation html"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-s", "--singleHtml", type=str, required=True, dest="singleHtmlFile", help="singleHtml file from sphinx-build")
parser.add_argument("-t", "--templateHtml", type=str, required=True, dest="templateHtmlFile", help="Template html file to insert docs into")

args = parser.parse_args()

soup = BeautifulSoup(open(args.singleHtmlFile), "lxml")

########################################
# Parse TOC
########################################

TOC = OrderedDict()
TOCMap = {}

# This marks the beginning of the TOC <div class="sphinxsidebarwrapper">

tocTag = soup.find("div", class_="sphinxsidebarwrapper")

# TOC entries look like this <li class="toctree-l1">

tocLevels = []

prevLevel = 1
prevName = ""

for li in tocTag.find_all("li"):
    
    name = ""
    
    for child in li.children:
        if child.name == "a":
            name = child['href']
            name = re.sub("index.html#","",name)
    
    if (prevLevel < int(li['class'][0][-1])):
        tocLevels.append(prevName)
        set(TOC, tocLevels, name)
        prevLevel = int(li['class'][0][-1])
        
    elif (prevLevel > int(li['class'][0][-1])):
        diff = prevLevel - int(li['class'][0][-1])
        for x in range(0, diff):
            tocLevels.pop()
        set(TOC, tocLevels, name)
        prevLevel = int(li['class'][0][-1])
        
    else :
        set(TOC, tocLevels, name)
        
    TOCMap[name] = int(li['class'][0][-1])
    
    prevName = name

########################################
# PRINT TOC
########################################

walkthroughTOC(soup, TOC, TOCMap)

