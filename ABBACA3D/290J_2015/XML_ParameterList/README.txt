This is a "stripped down" version of tinyxml 2.6.2 obtained from

http://www.grinninglizard.com/tinyxml/

and contains only those components of the distribution required to create the XML_ParameterListArrray
class.

The original source has been modified -- in tinyxml.h the lines 

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

were added to the header to always force the use of the STL string class.


May 29, 2011
