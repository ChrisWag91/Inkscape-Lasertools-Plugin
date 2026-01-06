# Inkscape-Lasertools-Plugin
A simple plugin to convert Inkscape vector graphics to Gcode for DIY laser engravers

## Version 1.x
The version found in the folder v1.0 is intended to be used with v1.x of Inkscape.
This Version is work in progress. Please feel free to test it, but don't rely on it in production.

## Version 0.9x
The version found in the folder v0.92 is intended to be used with v0.9x of Inkscape.
(Tested with Inkscape 0.92)

This project is based on the Gcodetools extension.

*******************************************************************************************************************************

## How to use Lasertools

**Watch this Video for a tutorial:**

<a href="https://youtu.be/F2U58onWOFM
" target="_blank"><img src="http://img.youtube.com/vi/F2U58onWOFM/0.jpg"
alt="Detailed Youtube Video" width="240" height="180" border="10" /></a>


**Watch this Video for some Backgrund:**

<a href="http://www.youtube.com/watch?feature=player_embedded&v=NhUvRJsa4D0
" target="_blank"><img src="http://img.youtube.com/vi/NhUvRJsa4D0/0.jpg"
alt="Detailed Youtube Video" width="240" height="180" border="10" /></a>


************************************************************************************
**step by step instructions:**

If you work from a picture (.png, .jpeg ...):

1. import picture to Inkscape
2. set the size of the page to the size you want your engraving to be [under File/Document Properties/Page]
3. **IMPORTANT:** set *Display Units* and *Units* to *mm* and set *Scale x* to 1 [under File/Document Properties/Page]
4. convert your picture to vector graphics [Path/Trace Bitmap]
5. remove the bitmap/png/jpeg... so that only the path is left
6. open Lasertools and set your parameters [under Extensions/Lasertools]
7. click Apply
8. **IMPORTANT:** the calculation takes some time, so be patient
9. finished

If you work from a svg file:

1. ungroup all objects until there are no more groups left
2. convert all objects to paths (especially text) using [Path/Object to Path] or [Path/Stroke to Path]
3. set the size of the page to the size you want your engraving to be [under File/Document Properties/Page]
4. **IMPORTANT:** set *Display Units* and *Units* to *mm* and set *Scale x* to 1 [under File/Document Properties/Page]
5. open Lasertools and set your parameters [under Extensions/Lasertools]
6. click Apply
7. **IMPORTANT:** the calculation takes some time, so be patient
8. finished

**Hint:**
If you end up with very complex paths (paths with thousands of points) you can simplify the path before exporting it to Gcode.
Use the function Path/Simplify.
The amount of simplification can be set under Edit/Preferences/Behavior.
A value of 0,0002 should work for most applications.

As of now, the first method seems to be more reliable.

Make sure the specified directory exists.
Calculating the Gcode can take quite some time on larger engravings. I would advice to test with a small engraving.

*******************************************************************************************************************************

## Installation
Copy the .py and .inx file from /v0.92 for Inkscape 0.9x or /v1.0 for Inkscape 1.x into the Inkscape extentions folder.

- under linux:    /usr/share/inkscape/extensions 
- under windows   C:\Program Files\Inkscape\share\extensions 

## Dependencies for v0.9x
- python2-lxml
- python2-numpy

## Dependencies for v1.0
- python-lxml
- python-numpy

**Hint:**

If you get an error on execution similar to: 
```console
"ModuleNotFoundError: No module named 'lxml'"
```
make sure to install the dependencies mentioned aboth.

On Linux for v1.0: 
```console
sudo apt-get install python-lxml
```
or
```console
sudo pacman -S python-lxml
```

*******************************************************************************************************************************

## Warning
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY Without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

License text: 
https://github.com/nevir/readable-licenses/blob/master/markdown/GPLv2-LICENSE.md

*******************************************************************************************************************************

## Open Bugs
- Fix infill glitches which appear on some paths
    - Workaround: make Path/Difference or Path Union operation which does not change the shape. That will fix the glitch.
    - Workaround: export as png, reimport and convert to paths again
- Scaling factor needs to be Set to 1
- Displayed units needs to be set to mm

## Features to be implemented / tested
- Automatically convert text to paths 
- Make a help section in the UI with the step by step instructions
- Performance improvements

*******************************************************************************************************************************

## Help appreciated
If you find bugs or want to implement a feature, please feel free to commit a change to the project.


