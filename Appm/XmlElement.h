#pragma once

#include <string>
#include <vector>
#include <ostream>

class XmlElement
{
public:
	XmlElement();
	XmlElement(const std::string & startTag, const std::string & endTag = "", const std::string & body = "");
	~XmlElement();

	void addChild(XmlElement * child);

	friend std::ostream & operator<<(std::ostream & os, const XmlElement & obj);

private:
	std::string startTag;
	std::string endTag;
	std::string body;
	std::vector<XmlElement*> children;
};

