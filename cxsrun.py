# ------------------------------------------------------------------------
#                           pyCXS 
#
#     GIMAP-CONICET Wind Turbine Cross Sectional Analysis Tool
# ------------------------------------------------------------------------
 
#  License
#  This file is part of pyCXS.
#
#  pyCXS is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  pyCXS is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pyCXS.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------


from code import cxs

filename = 'mh104'

xs = cxs.CrossSection(filename)

# Post-processing
xs.plotxs('stif', 6)
xs.plotxs('lyrs')
xs.plotxs('mlyr')
plt.axis('equal')
plt.show()
