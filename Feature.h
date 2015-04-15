#include <string>
//=============================================================================
using namespace std;
//=============================================================================
class Feature
{
private:
	string name;
	float value, min_value, max_value, incr_value;
public:
	Feature();
	Feature(string _name, float _value, float _min_value,
		float _max_value, float _incr_value);
	~Feature();

	string getName();
	void setName(string _name);
	float getValue();
	void setValue(float _value);
	float getMinValue();
	void setMinValue(float _min_value);
	float getMaxValue();
	void setMaxValue(float _max_value);
	float getIncrValue();
	void setIncrValue(float _incr_value);

	bool increment();
	bool decrement();
};