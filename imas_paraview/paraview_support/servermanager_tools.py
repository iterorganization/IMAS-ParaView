"""Helper functions for generating paraview ServerManager XMLs for the plugins."""

import functools
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
        entries: Key/value pairs to use as entries in the domain
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


def stringlistdomain(property_name, **kwargs):
    """Convenience decorator to add a StringListDomain to a property.

    Args:
        property_name: Name of the StringList property with allowed values.
        kwargs: Additional keyword arguments to add as xml attributes to the domain.
    """
    args = " ".join(f'{key}="{value}"' for key, value in kwargs.items())
    return smdomain.xml(
        f"<StringListDomain {args}><RequiredProperties>"
        f'<Property name="{property_name}" function="StringInfo" />'
        "</RequiredProperties></StringListDomain>"
    )


def arrayselectionstringvector(property_name, attribute_name):
    """Convenience decorator to add a stringvector for an ArraySelectionDomain.
    Paraview requires the following functions to be defined with the corresponding
    `attribute_name`:

        - GetNumberOf{attribute_name}Arrays(self)
        - Get{attribute_name}ArrayName(self, idx)
        - Get{attribute_name}ArrayStatus(self, *args)

    Args:
        property_name: Name of the stringvector property. Must match the
            `property_name` in arrayselectiondomain
        attribute_name: Name to use in the auxiliary Paraview functions
    """
    return smproperty.xml(
        f"""
        <StringVectorProperty
        information_only="1" name="{property_name}">
            <ArraySelectionInformationHelper attribute_name="{attribute_name}" />
        </StringVectorProperty>
        """
    )


def arrayselectiondomain(property_name, **kwargs):
    """Convenience decorator to add an ArraySelectionDomain to a property."""

    def decorator(func):
        args = " ".join(f'{key}="{value}"' for key, value in kwargs.items())

        xml = f"""
        <StringVectorProperty
            {args}
            command="{func.__name__}"
            number_of_elements="0"
            number_of_elements_per_command="2"
            panel_visibility="default"
            element_types="2 0"
            repeat_command="1">
            <ArraySelectionDomain name="array_list">
                <RequiredProperties>
                    <Property function="ArrayList" name="{property_name}" />
                </RequiredProperties>
            </ArraySelectionDomain>
            <Documentation>{func.__doc__}</Documentation>
            </StringVectorProperty>
        """

        return smproperty.xml(xml)(func)

    return decorator


def checkbox(**kwargs):
    """Convenience decorator for creating a simple boolean checkbox."""

    def decorator(func):
        args = " ".join(f'{key}="{value}"' for key, value in kwargs.items())

        xml = f"""
        <IntVectorProperty
            {args}
            command="{func.__name__}"
            number_of_elements="1">
            <BooleanDomain name="bool" />
            <Documentation>{func.__doc__}</Documentation>
        </IntVectorProperty>
        """

        return smproperty.xml(xml)(func)

    return decorator


def add_docstring(func):
    """Convenience decorator to add a Documentation XML node filled with the docstring
    of the property.
    """
    # Abuse smdomain.xml, which inserts the XML element in the correct location...
    doc = func.__doc__
    if doc:
        return smdomain.xml(f"<Documentation>{doc}</Documentation>")(func)
    return func


def _propertywrapper(property, *args, **kwargs):
    """Helper function to chain smproperty.<property> and add_docstring."""

    def decorator(func):
        return property(*args, **kwargs)(add_docstring(func))

    return decorator


intvector = functools.partial(_propertywrapper, smproperty.intvector)
"""Convenience decorator that wraps `smproperty.intvector` and `add_docstring`."""
stringvector = functools.partial(_propertywrapper, smproperty.stringvector)
"""Convenience decorator that wraps `smproperty.stringvector` and `add_docstring`."""
doublevector = functools.partial(_propertywrapper, smproperty.doublevector)
"""Convenience decorator that wraps `smproperty.doublevector` and `add_docstring`."""
