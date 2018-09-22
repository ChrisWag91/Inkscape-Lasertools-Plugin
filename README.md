# Inkscape-Lasertools-Plugin
A simple plugin to convert Inkscape vector graphics to Gcode for DIY laser engravers

## Version 0.91
Tested with Inkscape 0.92
This project is based on the Gcodetools extension.
*******************************************************************************************************************************

## How to use Lasertools
**Watch this Video for a quick introduction:**

<a href="http://www.youtube.com/watch?feature=player_embedded&v=NhUvRJsa4D0
" target="_blank"><img src="http://img.youtube.com/vi/NhUvRJsa4D0/0.jpg"
alt="Detailed Youtube Video" width="240" height="180" border="10" /></a>


************************************************************************************
**step by step instructions:**

1. import picture to Inkscape
2. set the size of the page to the size you want your engraving to be [under File/Document Properties/Page]
3. **IMPORTANT:** set *Display Units* and *Units* to *mm* and set *Scale x* to 1 [under File/Document Properties/Page]
4. convert your picture to vector graphics [Path/Trace Bitmap]
5. remove the bitmap/png/jpeg... so that only the path is left
6. open Lasertools and set your parameters [under Extensions/Lasertools]
7. click Apply
8. **IMPORTANT:** the calculation takes some time, so be patient
9. finished

Make sure the specified directory exists.
Only use the Live Preview for smaller engravings.
Calculating the Gcode can take quite some time on larger engravings. I would advice to test with a small engraving.

*******************************************************************************************************************************

## Installation
Copy the .py and .inx file to the Inkscape Extentions folder.

under linux:    /usr/share/inkscape/extensions
under windows   C:\Program Files\Inkscape\share\extensions 

*******************************************************************************************************************************

## Warning
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY Without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

License text: 
https://github.com/nevir/readable-licenses/blob/master/markdown/GPLv2-LICENSE.md

*******************************************************************************************************************************

## Open Bugs
- Scaling factor needs to be Set to 1
- displayed units needs to be set to mm
- sometimes small areas of an complicated path object get marked inversely

## Features to be implemented / tested
- Make Laser On/Off command configurable
- make offsets work
- make a help section in the ui with the step by step instructions
- performance improvements / multithreading
- implement checks for Scaling factor and display unit
- make text objects work

*******************************************************************************************************************************

## Help appreciated
If you find bugs or want to implement a feature, please feel free to commit a change to the project.
