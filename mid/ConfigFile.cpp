using namespace std;
#include <stdio.h>
#include <mxml.h>
#include<iostream>
#include<list>
#include<vector>
#include"ConfigFile.h"

/*
 * Customized xml reader for reading of user specified parameters.  
 * Expects files of the form:
 *
 * <?xml version="1.0" encoding="utf-8"?>
 * <Configuration>
 *   <ParametersGroup name="Movie Parameters">
 *     <IntegerParameter name="length" value="16384"/> 
 *     <IntegerParameter name="cx" value="4"/>
 *     <StringParameter name="file" value="~/data/file.dat"/> 
 *   </ParametersGroup>
 *   <ParametersGroup name="Fitting Parameters">
 *     <DoubleParameter name="tolerance" value="5e-5"/> 
 *   </ParametersGroup>
 * </Configuration>
 *
 *  See the test function at the end of the source for an example.
 *  Written by Matthew Grivich, The Salk Institute, October 2009.
 */

Parameter::Parameter(mxml_node_t *node)   
{
	this->node = node;
}

Parameter::Parameter(const Parameter &copyin)   
{       
   this->node = copyin.node;
}

const char *Parameter::getName() 
{
	return mxmlElementGetAttr(node, "name");
}

IntegerParameter::IntegerParameter(mxml_node_t *node) : Parameter(node) 
{
  mxml_node_t *currentNode = node;
  value = atoi(mxmlElementGetAttr(currentNode, "value")); 
}

IntegerParameter::IntegerParameter(const IntegerParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

int IntegerParameter::getValue() 
{
	return value;
}


IntegerListParameter::IntegerListParameter(mxml_node_t *node) : Parameter(node) 
{

  mxml_node_t *currentNode = node;
  char * buffer = (char *) malloc(strlen(mxmlElementGetAttr(currentNode, "value")) + 1);
   
  strcpy(buffer, mxmlElementGetAttr(currentNode, "value"));
  char * context = NULL;
  char * word = NULL;
  std::vector<string*> stringValue;
  for (word = strtok_r(buffer, ",", &context); word; word = strtok_r(NULL, ",", &context)) { 
	stringValue.push_back(new string(word)); 
  } 
  length = stringValue.size();

  value = (int *) malloc((stringValue.size()+1) * sizeof(int));
  value[0] = -1;
  for(int i=1; i<=stringValue.size(); i++) {
    value[i] = atoi(stringValue[i-1]->c_str());
  }
  
  while(!stringValue.empty()) {
    delete stringValue.back();
    stringValue.pop_back();
  }
   
}

IntegerListParameter::IntegerListParameter(const IntegerListParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

int* IntegerListParameter::getValue() 
{
	return value;
}

int IntegerListParameter::getIntegerCount() 
{
	return length;
}

IntegerListParameter::~IntegerListParameter() 
{
  free(value);
}

EnumeratorParameter::EnumeratorParameter(mxml_node_t *node) : Parameter(node) 
{
  mxml_node_t *currentNode = node;
  value = atoi(mxmlElementGetAttr(currentNode, "value")); 
}

EnumeratorParameter::EnumeratorParameter(const EnumeratorParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

int EnumeratorParameter::getValue() 
{
	return value;
}


DoubleParameter::DoubleParameter(mxml_node_t *node) : Parameter(node) 
{
  mxml_node_t *currentNode = node;
  value = atof(mxmlElementGetAttr(currentNode, "value"));
   
}

DoubleParameter::DoubleParameter(const DoubleParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

double DoubleParameter::getValue() 
{
	return value;
}



StringParameter::StringParameter(mxml_node_t *node) : Parameter(node) 
{
  mxml_node_t *currentNode = node;
  value = (char *) malloc(strlen(mxmlElementGetAttr(currentNode, "value")) + 1);
  strcpy(value, mxmlElementGetAttr(currentNode, "value"));
}

StringParameter::StringParameter(const StringParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

char * StringParameter::getValue() 
{
	return value;
}

StringParameter::~StringParameter() 
{
  free(value);
}


StringListParameter::StringListParameter(mxml_node_t *node) : Parameter(node) 
{
  mxml_node_t *currentNode = node;
  char * buffer = (char *) malloc(strlen(mxmlElementGetAttr(currentNode, "value")) + 1);
   strcpy(buffer, mxmlElementGetAttr(currentNode, "value"));
  char * context = NULL;
  char * word = NULL;
  for (word = strtok_r(buffer, ",", &context); word; word = strtok_r(NULL, ",", &context)) { 
	value.push_back(new string(word));
  } 
}

StringListParameter::StringListParameter(const StringListParameter &copyin) : Parameter(copyin) 
{ 
  this->value = copyin.value;
}

vector<string*> StringListParameter::getValue() 
{
	return value;
}

const char * StringListParameter::getCString(int i) {
  //does not make a copy.
  // -1 to match numerical recipes convention
	return value[i-1]->c_str(); 
}

int StringListParameter::getStringCount() 
{
	return value.size();
}

StringListParameter::~StringListParameter() 
{
  while(!value.empty()) {
    delete value.back();
    value.pop_back();
  }
}



ParameterGroup::ParameterGroup(mxml_node_t *node)  
{
  this->node = node;
  mxml_node_t *currentNode = node;
 
  //advance to the first parameter
  currentNode = mxmlWalkNext(node, node, MXML_DESCEND);
  currentNode = mxmlWalkNext(currentNode, node, MXML_DESCEND);

  
  while(currentNode != NULL) {
   
// This will print all the data in the file.   
// printf("Parameter Type: %s\n", currentNode->value.element.name);
// printf("Parameter Name: %s\n", mxmlElementGetAttr(currentNode, "name")); 
// printf("Parameter Value: %s\n\n", mxmlElementGetAttr(currentNode, "value"));
    
    char * paramType = currentNode->value.element.name;
	//skip if comment
    if(strncmp(paramType, "!--", 3) == 0);
    else if(strcmp(paramType, "IntegerParameter")==0)
       parameters.push_back(new IntegerParameter(currentNode));   
    else if(strcmp(paramType, "DoubleParameter")==0)
      parameters.push_back(new DoubleParameter(currentNode));
    else if(strcmp(paramType, "StringParameter")==0)
      parameters.push_back(new StringParameter(currentNode));
	else if(strcmp(paramType, "EnumeratorParameter")==0)
	  parameters.push_back(new EnumeratorParameter(currentNode));
    else if(strcmp(paramType, "StringListParameter")==0)
	  parameters.push_back(new StringListParameter(currentNode));
	else if(strcmp(paramType, "IntegerListParameter")==0)
	  parameters.push_back(new IntegerListParameter(currentNode));
	else
      cout << "Parameter Type \"" << paramType << "\" is not recognized.\n";
     
    
    //advance to the next parameter
    currentNode = mxmlWalkNext(currentNode, node, MXML_DESCEND);
    currentNode = mxmlWalkNext(currentNode, node, MXML_DESCEND);    
        
  }

}

ParameterGroup::ParameterGroup(const ParameterGroup &copyin)  
{       

   this->node = copyin.node;
   this->parameters = copyin.parameters;

}


const char *ParameterGroup::getName() 
{
  return mxmlElementGetAttr(node, "name");
}

Parameter *ParameterGroup::getParameter(const char *parameterName) 
{
  list<Parameter*>::iterator i;
  
  for(i=parameters.begin(); i!=parameters.end(); ++i) {
    Parameter *p = *i;
    if(strcmp(p->getName(), parameterName)==0) return p;
  }
  cout << "Parameter \"" << parameterName << "\" not found in parameter group \"" << this->getName() << ".\"\n";
  return NULL;
  
}

int ParameterGroup::getIntValue(const char *paramName) 
{
   IntegerParameter *ip = dynamic_cast<IntegerParameter*> (this->getParameter(paramName));
  if(ip == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type integer.\n";
  	return 0;
  }
  
  return ip->getValue();
}

int* ParameterGroup::getIntegerListValue(const char *paramName, int *count) 
{
  IntegerListParameter *dp = dynamic_cast<IntegerListParameter*> (this->getParameter(paramName));
  if(dp == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type integer list.\n";
	*count = -1;
  	return NULL;
  }
  *count = dp->getIntegerCount();
  return dp->getValue();
}

int ParameterGroup::getEnumValue(const char *paramName) 
{
   EnumeratorParameter *ip = dynamic_cast<EnumeratorParameter*> (this->getParameter(paramName));
  if(ip == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type enumerator.\n";
  	return 0;
  }
  
  return ip->getValue();
}

double ParameterGroup::getDoubleValue(const char *paramName) 
{
  DoubleParameter *dp = dynamic_cast<DoubleParameter*> (this->getParameter(paramName));
  if(dp == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type double.\n";
  	return 0.0;
  }
  
  return dp->getValue();
}

const char * ParameterGroup::getStringValue(const char *paramName) 
{
  StringParameter *dp = dynamic_cast<StringParameter*> (this->getParameter(paramName));
  if(dp == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type string.\n";
  	return "";
  }
  
  return dp->getValue();
}

/*
 * returns an array of strings.  Count is set to the number of strings. 
 */

char ** ParameterGroup::getStringListValue(const char *paramName, int *count) 
{
  StringListParameter *dp = dynamic_cast<StringListParameter*> (this->getParameter(paramName));
  if(dp == NULL) {
  	cout << "Parameter \"" << paramName << "\" is not of type string list.\n";
	*count = -1;
  	return NULL;
  }
  *count = dp->getStringCount();
  char ** strings = (char **) malloc((*count+1) * sizeof(char**));
  strings[0] = NULL;
  for(int i=1; i<=*count; i++) {
    const char *from = dp->getCString(i);
    char *to = (char *) malloc((strlen(from)+1));
	strcpy(to, from);
	strings[i] = to;
  }
  return strings;
}



ParameterGroup::~ParameterGroup() 
{
  while(!parameters.empty()) {
    delete parameters.back();
    parameters.pop_back();
  }
}



ConfigFile::ConfigFile(const char *fileName) 
{
  fp = fopen(fileName, "r");
  if(fp == 0) {
    cout << fileName << " not found.\n";
  }
  tree = mxmlLoadFile(NULL, fp, MXML_TEXT_CALLBACK);
  
  
  for (mxml_node_t *node = mxmlFindElement(tree, tree, "ParametersGroup", NULL, NULL, MXML_DESCEND); node!=NULL;
    node = mxmlFindElement(node, tree, "ParametersGroup", NULL, NULL, MXML_DESCEND) ) { 
    
    parameterGroups.push_back(new ParameterGroup(node));
  }
  
  fclose(fp);
   
}

ConfigFile::ConfigFile(const ConfigFile &copyin) 
{ //copy constructor to handle pass by value
  this->parameterGroups = copyin.parameterGroups;
  this->fp = copyin.fp;
  this->tree = copyin.tree;
}


ConfigFile::~ConfigFile() 
{
  mxmlDelete(tree); //Deletes entire tree.  This is a C style structure.
  
  while(!parameterGroups.empty()) {
    delete parameterGroups.back();
    parameterGroups.pop_back();
  }
}

ParameterGroup * ConfigFile::getParameterGroup(const char *groupName) 
{
  list<ParameterGroup*>::iterator i;
  for(i=parameterGroups.begin(); i!=parameterGroups.end(); ++i) {
    if(strcmp((**i).getName(), groupName) == 0) {
    	return *i;
    }
   
  }  
  cout << "Parameter Group \"" << groupName << "\" was not found.\n";
  return NULL;
}

int XMLTestFunction() 
{
  ConfigFile *cf = new ConfigFile("params.xml");
  
  ParameterGroup *pg = cf->getParameterGroup("Movie Parameters"); 
  printf("length: %d\n", pg->getIntValue("length"));
  printf("file: %s\n", pg->getStringValue("file"));
 
  pg = cf->getParameterGroup("Fitting Parameters"); 
  printf("tolerance: %f\n", pg->getDoubleValue("tolerance"));
  
  delete cf; 
  
  return 0;
}
