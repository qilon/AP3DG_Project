#include <string>
//=============================================================================
using namespace std;
//=============================================================================
struct FeatureConfig
{
	string name;
	float init_value, min_value, max_value, incr_value;

	FeatureConfig()
	{
		name = "";
		init_value = 0.f;
		min_value = numeric_limits<float>::min();
		max_value = numeric_limits<float>::max();
		incr_value = 1.f;
	};
	FeatureConfig(string _name, float _init_value, float _min_value, float _max_value,
		float _incr_value)
	{
		name = _name;
		init_value = _init_value;
		min_value = _min_value;
		max_value = _max_value;
		incr_value = _incr_value;
	};
};