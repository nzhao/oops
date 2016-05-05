#include "include/misc/xmlreader.h"

ConfigXML::ConfigXML(string filename)
{
    xml_document<> doc;
    xml_node<> *   root_node;

    ifstream theFile ( filename.c_str() );
    vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
    buffer.push_back('\0');

    doc.parse<0>(&buffer[0]);
    root_node = doc.first_node();

    for (xml_node<> * section = root_node->first_node(); section; section = section->next_sibling())
    {
        int count = 0;
        string section_name = section->name();
        for(xml_node<> * item = section->first_node(); item; item = item->next_sibling())
        {
            string item_name =  item->name();
            string item_type = item->first_attribute("type")->value();
            string item_value = item->value();
            _parameters[make_pair(section_name, item_name)] = make_pair(item_type, item_value);
            count++;
        }
    }

}

void ConfigXML::printParameters() const
{
    PARA_MAP::const_iterator pos;
    cout << endl;
    cout << "####################################### parameters #######################################" << endl;
    for(pos = _parameters.begin(); pos != _parameters.end(); ++pos)
    {
        cout <<  std::right << setw(12) << (pos->first).first << " :: " << std::left << setw(25) <<  (pos->first).second;
        cout <<  std::right << setw(12) << (pos->second).first << " = " ;
        cout <<  std::left << setw(35) << (pos->second).second << endl; 
    }
    cout << "##########################################################################################" << endl;
    cout << endl;
}

int ConfigXML::getIntParameter(string section_name, string para_name) const
{
    pair<string, string> name = make_pair(section_name, para_name);
    if( !strcmp(_parameters[name].first.c_str(), "int") )
        return atoi(_parameters[name].second.c_str());
    else
    {
        cout  << "Wrong parameter type: " << name.first << ":" << name.second << " is not INT type" << endl;;
        return DEFAULT_INT_PARAMETER;
    }
}

double ConfigXML::getDoubleParameter(string section_name, string para_name) const
{
    pair<string, string> name = make_pair(section_name, para_name);
    if( !strcmp(_parameters[name].first.c_str(), "double") )
        return atof(_parameters[name].second.c_str());
    else
    {
        cout  << "Wrong parameter type: " << name.first << ":" << name.second << " is not DOUBLE type" << endl;;
        return DEFAULT_DOUBLE_PARAMETER;
    }
}

string ConfigXML::getStringParameter(string section_name, string para_name) const
{
    pair<string, string> name = make_pair(section_name, para_name);
    if( !strcmp(_parameters[name].first.c_str(), "char") )
        return _parameters[name].second;
    else
    {
        cout  << "Wrong parameter type: " << name.first << ":" << name.second << " is not CHAR type" << endl;;
        return DEFAULT_STRING_PARAMETER;
    }
}

vec ConfigXML::getVectorParameter(string section_name, string para_name) const
{
    pair<string, string> name = make_pair(section_name, para_name);
    if( !strcmp(_parameters[name].first.c_str(), "vector") )
    {
        string vString =  _parameters[name].second;
        return vec(vString);
    }
    else
    {
        cout  << "Wrong parameter type: " << name.first << ":" << name.second << " is not CHAR type" << endl;;
        return DEFAULT_STRING_PARAMETER;
    }
}
