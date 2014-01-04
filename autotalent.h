#ifndef __AUTOTALENT__
#define __AUTOTALENT__

#include "IPlug_include_in_plug_hdr.h"

class autotalent : public IPlug
{
public:
  autotalent(IPlugInstanceInfo instanceInfo);
  ~autotalent();

  void Reset();
  void OnParamChange(int paramIdx);
  void ProcessDoubleReplacing(double** inputs, double** outputs, int nFrames);

private:
  double mGain;
};

#endif
