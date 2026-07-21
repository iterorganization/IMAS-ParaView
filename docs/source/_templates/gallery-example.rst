.. _{{ entry.name }}:

{{ entry.title }}
{{ '=' * entry.title|length }}

.. figure:: /gallery/examples/{{ entry.dir_name }}/{{ entry.image_name }}
   :alt: {{ entry.title }}
   :align: center

**Author:** {{ entry.author }}

{{ entry.description }}

**Data URI(s):**

{% for uri in entry.uris %}
.. code-block:: text

   {{ uri }}
{% endfor %}


{% if entry.link %}
**Link:** `View resource <{{ entry.link }}>`_
{% endif %}

{% if entry.version %}
**IMAS-ParaView version:** `{{ entry.version }} <https://github.com/iterorganization/IMAS-ParaView/releases/tag/{{ entry.version }}>`_
{% endif %}

{% if entry.state_file_name %}
:download:`Download ParaView State File <{{ entry.state_file_name }}>`
{% endif %}
