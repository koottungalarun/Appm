#include "XmlElement.h"



XmlElement::XmlElement()
{
}

XmlElement::XmlElement(const std::string & startTag, const std::string & endTag, const std::string & body)
{
	this->startTag = startTag;
	this->endTag = endTag;
	this->body = body;
}


XmlElement::~XmlElement()
{
	if (children.size() > 0) {
		for (auto child : children) {
			delete child;
		}
	}
	children.resize(0);
}

void XmlElement::addChild(XmlElement * child)
{
	children.push_back(child);
}

std::ostream & operator<<(std::ostream & os, const XmlElement & obj)
{
	os << obj.startTag << std::endl;
	if (obj.body.size() > 0) {
		os << obj.body << std::endl;
	}
	for (auto child : obj.children) {
		os << *child << std::endl;
	}
	os << obj.endTag;
	return os;
}
