# -*- coding: utf-8 -*-
#
# Copyright (c) 2020 Martin Owens <doctormo@gmail.com>
#                    Thomas Holder <thomas.holder@schrodinger.com>
#                    Sergei Izmailov <sergei.a.izmailov@gmail.com>
#                    Windell Oskay <windell@oskay.net>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# pylint: disable=attribute-defined-outside-init
#
"""
Provide a way to load lxml attributes with an svg API on top.
"""

import random
from collections import OrderedDict
from lxml import etree

from ..units import discover_unit, convert_unit, render_unit
from ..transforms import BoundingBox
from ..styles import StyleSheets

from ._base import BaseElement
from ._meta import NamedView, Defs, StyleElement, Metadata

if False: # pylint: disable=using-constant-test
    import typing # pylint: disable=unused-import


class SvgDocumentElement(BaseElement): # pylint: disable=too-many-public-methods
    """Provide access to the document level svg functionality"""
    tag_name = 'svg'

    def _init(self):
        self.current_layer = None
        self.view_center = (0.0, 0.0)
        self.selected = OrderedDict()
        self.ids = {}

    def tostring(self):
        """Convert document to string"""
        return etree.tostring(etree.ElementTree(self))

    def get_ids(self):
        """Returns a set of unique document ids"""
        if not self.ids:
            self.ids = set(self.xpath('//@id'))
        return self.ids

    def get_unique_id(self, prefix, size=4):
        """Generate a new id from an existing old_id"""
        ids = self.get_ids()
        new_id = None
        _from = 10 ** size - 1
        _to = 10 ** size
        while new_id is None or new_id in ids:
            # Do not use randint because py2/3 incompatibility
            new_id = prefix + str(int(random.random() * _from - _to) + _to)
        self.ids.add(new_id)
        return new_id

    def set_selected(self, *ids):
        """
        Sets the currently selected elements to these ids.

        Arguments a list of element ids, element objects or
            a single xpath expression starting with "//".

        All element objects must have an id to be correctly set.

        >>> svg.set_selected("rect123", "path456", "text789")
        >>> svg.set_selected(elem1, elem2, elem3)
        >>> svg.set_selected("//rect")
        """
        self.selected = OrderedDict()

        # Allow selecting of xpath elements directly
        if len(ids) == 1 and isinstance(ids[0], str) and ids[0].startswith('//'):
            ids = self.xpath(ids[0])

        for elem_id in ids:
            if isinstance(elem_id, BaseElement):
                # Selection is a list of nodes to select
                self.selected[elem_id.get('id')] = elem_id
                continue
            # Selection is a text element id, find it (or them).
            for node in self.xpath('//*[@id="{}"]'.format(elem_id)):
                self.selected[elem_id] = node

    def get_z_selected(self):
        """Get the selected elements, but ordered by their apperence in the document"""
        sel = self.selected
        return OrderedDict((_id, sel[_id]) for _id in self.xpath('//@id') if _id in sel)

    def get_selected(self, *types):
        """Generator: Gets selected nodes which are the given element types"""
        for node in self.selected.values():
            if not types or isinstance(node, types):
                yield node

    def get_selected_or_all(self, *types):
        """Returns a generator of selected items: i.e. svg.get_selected(types)
             or all of this type of element i.e. svg.descendants(types)
        """
        if self.selected:
            for node in self.get_selected(*types):
                yield node # yield from when py3 only
        else:
            for node in self.descendants(*types):
                yield node # yield from when py3 only

    def get_selected_bbox(self):
        """
        Gets a :class:`inkex.transforms.BoundingBox` object for the selected items.

        Text objects have a bounding box without width or height that only
        reflects the coordinate of their anchor. If a text object is a part of
        the selection's boundary, the bounding box may be inaccurate.

        When no object is selected or when the object's location cannot be
        determined (e.g. empty group or layer), all coordinates will be None.
        """
        return sum([node.bounding_box() for node in self.selected.values()], None)

    def get_page_bbox(self):
        """Gets the page dimensions as a bbox"""
        return BoundingBox((0, float(self.width)), (0, float(self.height)))

    def get_first_selected(self, *types):
        """Returns the first item in the selected list, of the given types"""
        if self.selected:
            return list(self.get_selected(*types))[0]
        return None

    def get_current_layer(self):
        """Returns the currently selected layer"""
        layer = self.getElementById(self.namedview.current_layer, 'svg:g')
        if layer is None:
            return self
        return layer

    def getElement(self, xpath):  # pylint: disable=invalid-name
        """Gets a single element from the given xpath or returns None"""
        return self.findone(xpath)

    def getElementById(self, eid, elm='*'):  # pylint: disable=invalid-name
        """Get an element in this svg document by it's ID attribute"""
        if eid is not None:
            eid = eid.strip()[4:-1] if eid.startswith('url(') else eid
            eid = eid.lstrip('#')
        return self.getElement('//{}[@id="{}"]'.format(elm, eid))

    def getElementsByHref(self, eid): # pylint: disable=invalid-name
        """Get elements by their href xlink attribute"""
        return self.xpath('//*[@xlink:href="#{}"]'.format(eid))

    def getElementsByStyleUrl(self, eid, style=None): # pylint: disable=invalid-name
        """Get elements by a style attribute url"""
        url = "url(#{})".format(eid)
        if style is not None:
            url = style + ":" + url
        return self.xpath('//*[contains(@style,"{}")]'.format(url))

    @property
    def name(self):
        """Returns the Document Name"""
        return self.get('sodipodi:docname', '')

    @property
    def namedview(self):
        """Return the sp namedview meta information element"""
        return self.get_or_create('//sodipodi:namedview', NamedView, True)

    @property
    def metadata(self):
        """Return the svg metadata meta element container"""
        return self.get_or_create('//svg:metadata', Metadata, True)

    @property
    def defs(self):
        """Return the svg defs meta element container"""
        return self.get_or_create('//svg:defs', Defs, True)

    def get_viewbox(self):
        """Parse and return the document's viewBox attribute"""
        try:
            ret = [float(unit) for unit in self.get('viewBox', '0').split()]
        except ValueError:
            ret = ''
        if len(ret) != 4:
            return [0, 0, 0, 0]
        return ret

    @property
    def width(self):  # getDocumentWidth(self):
        """Fault tolerance for lazily defined SVG"""
        return self.unittouu(self.get('width')) or self.get_viewbox()[2]

    @property
    def height(self):  # getDocumentHeight(self):
        """Returns a string corresponding to the height of the document, as
        defined in the SVG file. If it is not defined, returns the height
        as defined by the viewBox attribute. If viewBox is not defined,
        returns the string '0'."""
        return self.unittouu(self.get('height')) or self.get_viewbox()[3]

    @property
    def scale(self):
        """Return the ratio between the page width and the viewBox width"""
        try:
            scale_x = float(self.width) / float(self.get_viewbox()[2])
            scale_y = float(self.height) / float(self.get_viewbox()[3])
            return max([scale_x, scale_y])
        except (ValueError, ZeroDivisionError):
            return 1.0

    @property
    def unit(self):
        """Returns the unit used for in the SVG document.
        In the case the SVG document lacks an attribute that explicitly
        defines what units are used for SVG coordinates, it tries to calculate
        the unit from the SVG width and viewBox attributes.
        Defaults to 'px' units."""
        viewbox = self.get_viewbox()
        if viewbox and set(viewbox) != {0}:
            return discover_unit(self.get('width'), viewbox[2], default='px')
        return 'px'  # Default is px

    def unittouu(self, value):
        """Convert a unit value into the document's units"""
        return convert_unit(value, self.unit)

    def uutounit(self, value, to_unit):
        """Convert from the document's units to the given unit"""
        return convert_unit(render_unit(value, self.unit), to_unit)

    def add_unit(self, value):
        """Add document unit when no unit is specified in the string """
        return render_unit(value, self.unit)

    @property
    def stylesheets(self):
        """Get all the stylesheets, bound together to one, (for reading)"""
        sheets = StyleSheets(self)
        for node in self.xpath('//svg:style'):
            sheets.append(node.stylesheet())
        return sheets

    @property
    def stylesheet(self):
        """Return the first stylesheet or create one if needed (for writing)"""
        for sheet in self.stylesheets:
            return sheet

        style_node = StyleElement()
        self.defs.append(style_node)
        return style_node.stylesheet()
