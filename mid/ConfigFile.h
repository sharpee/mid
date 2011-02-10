#include<mxml.h>
#include<iostream>
#include<string>
#include<list>
#include<vector>

/*
 * Header file for ConfigFile.cpp.  Look there for the comments.
 */

class Parameter 
{
   virtual void placeholder(){}; // This allows the class to be extended 
                                 // and for dynamic_cast to work.

   public:
      mxml_node_t *node;

      Parameter(mxml_node_t *node);
      Parameter(const Parameter &);
      ~Parameter(){};


  const char * getName();
};

class IntegerParameter : public Parameter  
{
  int value;

  public:
    IntegerParameter(mxml_node_t *node);
    IntegerParameter(const IntegerParameter &);
    ~IntegerParameter(){};
    
    int getValue();
};

class IntegerListParameter : public Parameter  
{
  int *value;
  int length;

  public:
    IntegerListParameter(mxml_node_t *node);
    IntegerListParameter(const IntegerListParameter &);
	int getInteger(int i);
	int getIntegerCount();
    ~IntegerListParameter();
    
    int* getValue();
};

class EnumeratorParameter : public Parameter  
{
  int value;

  public:
    EnumeratorParameter(mxml_node_t *node);
    EnumeratorParameter(const EnumeratorParameter &);
    ~EnumeratorParameter(){};
    
    int getValue();
};

class DoubleParameter : public Parameter  
{
  double value;

  public:
    DoubleParameter(mxml_node_t *node);
    DoubleParameter(const DoubleParameter &);
    ~DoubleParameter(){};
    
    double getValue();
};

class StringParameter : public Parameter  
{
  char * value;

  public:
    StringParameter(mxml_node_t *node);
    StringParameter(const StringParameter &);
    ~StringParameter();
    
    char * getValue();
};

class StringListParameter : public Parameter  
{
  std::vector<string*> value;

  public:
    StringListParameter(mxml_node_t *node);
    StringListParameter(const StringListParameter &);
	const char * getCString(int i);
	int getStringCount();
    ~StringListParameter();
    
    std::vector<string*> getValue();
};

class ParameterGroup 
{
   public:

    mxml_node_t *node;
    std::list<Parameter*> parameters;

    ParameterGroup(mxml_node_t *node);
    ParameterGroup(const ParameterGroup &);
    ~ParameterGroup();

    const char * getName();
    Parameter *getParameter(const char *name);
    int getIntValue(const char *paramName);
	int* getIntegerListValue(const char *paramName, int *count);
	int getEnumValue(const char *paramName);
    double getDoubleValue(const char *paramName);
    const char * getStringValue(const char *paramName);
	char ** getStringListValue(const char *paramName, int *count);
  
};

class ConfigFile  
{
  list<ParameterGroup*> parameterGroups;
  FILE *fp;
  mxml_node_t *tree;

  public:
    ConfigFile(const char *fileName);
    ConfigFile(const ConfigFile &);
    ~ConfigFile();
    ParameterGroup *getParameterGroup(const char *groupName);
   
};

