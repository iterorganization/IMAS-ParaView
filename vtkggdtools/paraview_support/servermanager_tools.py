"""Helper functions for generating paraview ServerManager XMLs for the plugins.
"""

from typing import List

from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty


def propertygroup(label: str, properties: List[str]):
    """Convenience decorator for creating a property group.

    Args:
        label: Label to use for the property group
        properties: Names of the properties to include in the group
    """
    prop_xml = "".join(f'<Property name="{prop}" />' for prop in properties)
    return smproperty.xml(f'<PropertyGroup label="{label}">{prop_xml}</PropertyGroup>')


def enumeration(name: str, entries: dict):
    """Convenience decorator for creating an EnumerationDomain.

    Args:
        name: Name of the EnumerationDomain
        kwargs: Key/value pairs to use as entries in the domain
    """
    vals = "".join(f'<Entry text="{k}" value="{v}" />' for k, v in entries.items())
    return smdomain.xml(f'<EnumerationDomain name="{name}">{vals}</EnumerationDomain>')


def genericdecorator(**kwargs):
    """Convenience decorator to add a PropertyWidgetDecorator hint.

    This can be used to hide or disable property widgets based on the value of another
    property.
    """
    args = " ".join(f'{key}="{value}"' for key, value in kwargs.items())
    return smhint.xml(f'<PropertyWidgetDecorator type="GenericDecorator" {args}/>')
