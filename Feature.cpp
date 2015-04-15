#include "Feature.h"
//=============================================================================
Feature::Feature()
{
	name = "";
	value = 0.f;
	min_value = 0.f;
	max_value = 0.f;
	incr_value = 0.f;
}
//=============================================================================
Feature::Feature(string _name, float _value, float _min_value,
	float _max_value, float _incr_value)
{
	name = _name;
	value = _value;
	min_value = _min_value;
	max_value = _max_value;
	incr_value = _incr_value;
}
//=============================================================================
Feature::~Feature()
{

}
//=============================================================================
string Feature::getName()
{
	return name;
}
//=============================================================================
void Feature::setName(string _name)
{
	name = _name;
}
//=============================================================================
float Feature::getValue()
{
	return value;
}
//=============================================================================
void Feature::setValue(float _value)
{
	value = _value;
}
//=============================================================================
float Feature::getMinValue()
{
	return min_value;
}
//=============================================================================
void Feature::setMinValue(float _min_value)
{
	min_value = _min_value;
}
//=============================================================================
float Feature::getMaxValue()
{
	return max_value;
}
//=============================================================================
void Feature::setMaxValue(float _max_value)
{
	max_value = _max_value;
}
//=============================================================================
float Feature::getIncrValue()
{
	return incr_value;
}
//=============================================================================
void Feature::setIncrValue(float _incr_value)
{
	incr_value = _incr_value;
}
//=============================================================================
bool Feature::increment()
{
	float new_value = value + incr_value;
	if (new_value > max_value)
	{
		return false;
	}
	else
	{
		value = new_value;
		return true;
	}
}
//=============================================================================
bool Feature::decrement()
{
	float new_value = value - incr_value;
	if (new_value < min_value)
	{
		return false;
	}
	else
	{
		value = new_value;
		return true;
	}
}
//=============================================================================