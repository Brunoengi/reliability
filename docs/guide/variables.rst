Variables Distribution 
======================

This section was created to help with the types of random variables and their parameters. This repository works with eight types of random variables, five of which are called usual variables or typical variables and three are called extreme value variables.

Typical Distribution
---------------------------------------------------

- `Normal Distribution <#id1>`_
- `Log-normal Distribution <#id2>`_
- `Uniform Distribution <#id3>`_
- `Beta Distribution <#id4>`_
- `Gamma Distribution <#id5>`_


Extreme Value Distribution
---------------------------------------------------
- `Gumbel Distribution <#id6>`_
- `Frechet Distribution <#id7>`_
- `Weibull Distribution <#id8>`_

How to create a variable set
--------------------------------------------------

Variables must be entered in the form of a list where each position in the list will be a dictionary. In the following example, we have a problem involving three variables.
  .. code-block:: bash

    xvar = [
      {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10},
      {'varname': 'Z', 'vardist': 'normal', 'varmean': 60.00, 'varcov': 0.05},
      {'varname': 'M', 'vardist': 'frechet', 'varmean': 1000.00, 'varcov': 0.30}
    ]

However, the common mandatory parameters of all probability distributions are:
``varname``,
``vardist``

The other dictionary keys depend on each distribution type. In this case, it is necessary to read the documentation of the distribution of interest.
  

Normal Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.


Log-normal Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.

Uniform Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.

Beta Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.

Gamma Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.


Gumbel Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.


Frechet Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.

Weibull Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas et convallis ligula, aliquam consectetur nunc. Aliquam at justo non sapien euismod lobortis. In lacinia fringilla semper. Sed sed ex sit amet felis aliquet congue quis vel turpis. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas vitae gravida dui, quis consectetur est. Donec sodales magna sed nisi viverra commodo. Nullam sed lectus euismod, blandit libero eu, ornare nulla. Nunc congue fermentum metus, sit amet elementum leo vestibulum ac. Duis vel congue elit. Integer porttitor tellus nec feugiat malesuada. Ut cursus sapien ac sapien suscipit, nec rutrum lacus semper.
