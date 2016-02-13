#include "include/misc/xmlreader.h"

ConfigXML::ConfigXML(string filename)
{
    //xml_document<> _doc;
    //xml_node<> * root_node;
    ifstream theFile (filename);
    vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
    buffer.push_back('\0');

    _doc.parse<0>(&buffer[0]);
    _root_node = _doc.first_node();

    for (xml_node<> * section = _root_node->first_node(); section; section = section->next_sibling())
    {
        int count = 0;
        string section_name = section->name();
        for(xml_node<> * item = section->first_node(); item; item = item->next_sibling())
        {
            string item_name =  item->name();
            string item_type = item->first_attribute("type")->value();
            string item_value = item->value();
            _parameters[item_name] = make_pair(item_type, item_value);
            count++;
        }
        cout << "Section: " << section_name << ", " << count << " parameters read." << endl;
    }

}

void ConfigXML::printParameters()
{
    PARA_MAP::iterator pos;
    for(pos = _parameters.begin(); pos != _parameters.end(); ++pos)
    {
        cout << "para: " << pos->first << "\t" ;
        cout << "type: " << (pos->second).first << "\t";
        cout << "value: " << (pos->second).second << endl; 
    }
}

int ConfigXML::getIntParameter(string name)
{
    if( !strcmp(_parameters[name].first.c_str(), "int") )
        return atoi(_parameters[name].second.c_str());
    else
    {
        cout  << "Wrong parameter type: " << name << " is not INT type" << endl;;
        return DEFAULT_INT_PARAMETER;
    }
}

double ConfigXML::getDoubleParameter(string name)
{
    if( !strcmp(_parameters[name].first.c_str(), "double") )
        return atof(_parameters[name].second.c_str());
    else
    {
        cout  << "Wrong parameter type: " << name << " is not DOUBLE type" << endl;;
        return DEFAULT_DOUBLE_PARAMETER;
    }
}

string ConfigXML::getStringParameter(string name)
{
    if( !strcmp(_parameters[name].first.c_str(), "char") )
        return _parameters[name].second;
    else
    {
        cout  << "Wrong parameter type: " << name << " is not CHAR type" << endl;;
        return DEFAULT_STRING_PARAMETER;
    }
}

