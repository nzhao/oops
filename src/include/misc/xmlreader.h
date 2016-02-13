#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "include/rapidxml-1.13/rapidxml.hpp"

using namespace rapidxml;
using namespace std;

typedef map<string, pair<string, string> > PARA_MAP;
const int     DEFAULT_INT_PARAMETER = 0;
const double  DEFAULT_DOUBLE_PARAMETER = 0.0;
const string  DEFAULT_STRING_PARAMETER = "None";

class ConfigXML
{
public:
    ConfigXML() {};
    ConfigXML(string filename);
    ~ConfigXML() {};

    void   printParameters();
    int    getIntParameter(string para_name);
    double getDoubleParameter(string para_name);
    string getStringParameter(string para_name);

protected:
private:
    xml_document<> _doc;
    xml_node<> *   _root_node;
    PARA_MAP       _parameters;
};

//int main(void)
//{
    //string filename = "config.xml";
    //xml_document<> doc;
    //xml_node<> * root_node;
    //Read the xml file into a vector
        //ifstream theFile (filename);
    //vector<char> buffer((istreambuf_iterator<char>(theFile)), istreambuf_iterator<char>());
    //buffer.push_back('\0');

    //// Parse the buffer using the xml file parsing library into doc 
    //doc.parse<0>(&buffer[0]);
    //// Find our root node
    //root_node = doc.first_node();

    //PARA_MAP parameters;
    //// Iterate over the brewerys
    //for (xml_node<> * section = root_node->first_node(); section; section = section->next_sibling())
    //{
        //int count = 0;
        //string section_name = section->name();
        //for(xml_node<> * item = section->first_node(); item; item = item->next_sibling())
        //{
            //string item_name =  item->name();
            //string item_type = item->first_attribute("type")->value();
            //string item_value = item->value();
            //parameters[item_name] = make_pair(item_type, item_value);
            //count++;
        //}
        //cout << "Section: " << section_name << ", " << count << " parameters read." << endl;
    //}

    //PARA_MAP::iterator pos;
    //for(pos = parameters.begin(); pos != parameters.end(); ++pos)
    //{
        //cout << "para: " << pos->first << "\t" ;
        //cout << "type: " << (pos->second).first << "\t";
        //cout << "value: " << (pos->second).second << endl; 
    //}
    //cout << parameters["method"].second << endl;
//}

