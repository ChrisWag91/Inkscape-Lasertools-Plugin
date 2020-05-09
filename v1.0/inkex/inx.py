# coding=utf-8
#
# Copyright (C) 2018 Martin Owens <doctormo@gmail.com>
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
"""
Parsing inx files for checking and generating.
"""

import os
from lxml import etree

try:
    from inspect import isclass
    from importlib import util
except ImportError:
    util = None  # type: ignore

from .base import InkscapeExtension
from .utils import Boolean

NSS = {
    'inx': 'http://www.inkscape.org/namespace/inkscape/extension',
}

class InxLookup(etree.CustomElementClassLookup):
    """Custom inx xml file lookup"""
    def lookup(self, node_type, document, namespace, name):  # pylint: disable=unused-argument
        if name == 'param':
            return ParamElement
        return None

INX_PARSER = etree.XMLParser()
INX_PARSER.set_element_class_lookup(InxLookup())

class InxFile(object):
    """Open an INX file and provide useful functions"""
    name = property(lambda self: self._text('inx:name'))
    ident = property(lambda self: self._text('inx:id'))
    slug = property(lambda self: self.ident.split('.')[-1].title().replace('_', ''))
    kind = property(lambda self: self.metadata['type'])

    def __init__(self, filename):
        self.filename = os.path.basename(filename)
        self.doc = etree.parse(filename, parser=INX_PARSER)
        self.root = self.doc.getroot()

    def __repr__(self):
        return "<inx '{0.filename}' '{0.name}'>".format(self)

    def xpath(self, xpath):
        """Namespace specific xpath searches"""
        return self.root.xpath(xpath, namespaces=NSS)

    def find_one(self, name):
        """Return the first element matching the given name"""
        for elem in self.xpath(name):
            return elem
        return None

    def _text(self, name, default=None):
        elem = self.find_one(name)
        if elem is not None and elem.text:
            return elem.text
        return default

    @property
    def script(self):
        """Returns information about the called script"""
        command = self.find_one('inx:script/inx:command')
        if command is None:
            return {}
        return {
            'interpreter': command.get('interpreter', None),
            'location': command.get('location', None),
            'script': command.text,
        }

    @property
    def extension_class(self):
        """Attempt to get the extension class"""
        script = self.script.get('script', None)
        if script is not None and util is not None:
            name = script[:-3].replace('/', '.')
            spec = util.spec_from_file_location(name, script)
            mod = util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            for value in mod.__dict__.values():
                if 'Base' not in name and isclass(value) and  value.__module__ == name \
                    and issubclass(value, InkscapeExtension):
                    return value
        return None

    @property
    def metadata(self):
        """Returns information about what type of extension this is"""
        effect = self.find_one('inx:effect')
        output = self.find_one('inx:output')
        data = {}
        if effect is not None:
            data['type'] = 'effect'
            data['preview'] = Boolean(effect.get('needs-live-preview', 'true'))
            data['objects'] = self._text('inx:effect/inx:object-type', 'all')
        elif self.find_one('inx:input') is not None:
            data['type'] = 'input'
            data['extension'] = self._text('inx:input/inx:extension')
            data['mimetype'] = self._text('inx:input/inx:mimetype')
            data['name'] = self._text('inx:input/inx:filetypename')
            data['tooltip'] = self._text('inx:input/inx:filetypetooltip')
        elif output is not None:
            data['type'] = 'output'
            data['dataloss'] = Boolean(self._text('inx:output/inx:dataloss', 'false'))
            data['extension'] = self._text('inx:output/inx:extension')
            data['mimetype'] = self._text('inx:output/inx:mimetype')
            data['name'] = self._text('inx:output/inx:filetypename')
            data['tooltip'] = self._text('inx:output/inx:filetypetooltip')
        return data

    @property
    def menu(self):
        """Return the menu this effect ends up in"""
        def _recurse_menu(parent):
            for child in parent.xpath('inx:submenu', namespaces=NSS):
                yield child.get('name')
                for subchild in _recurse_menu(child):
                    yield subchild
                break # Not more than one menu chain?
        menu = self.find_one('inx:effect/inx:effects-menu')
        return list(_recurse_menu(menu)) + [self.name]

    @property
    def params(self):
        """Get all params at all levels"""
        # Returns any params at any levels
        return list(self.xpath('//inx:param'))


class ParamElement(etree.ElementBase):
    """
    A param in an inx file.
    """
    name = property(lambda self: self.get('name'))
    param_type = property(lambda self: self.get('type', 'string'))

    @property
    def options(self):
        """Return a list of option values"""
        if self.param_type == 'notebook':
            return [option.get('name') for option in self.xpath('inx:page', namespaces=NSS)]
        return [option.get('value') for option in self.xpath('inx:option', namespaces=NSS)]

    def __repr__(self):
        return "<param name='{0.name}' type='{0.param_type}'>".format(self)
