#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

namespace
{
  float normAng(const float x)
  {
    float val = fmodf(x + F_PI, 2.0F*F_PI);
    if (val < 0.0F)
    {val += 2.0F*F_PI;}
    return val - F_PI;
  }
}

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{ const float len = L*0.5F*sqrt(2.0F);
  const float len_inv = 1.0F / len;
  const float forc_thrust = collThrustCmd;
  const float TX = momentCmd.x, TY = momentCmd.y, TZ = momentCmd.z;
  const float kVariable_inverse = 1.0F / kappa;  
  cmd.desiredThrustsN[0] = (1/(4.0)) * (forc_thrust + len_inv*TX + len_inv*TY - kVariable_inverse*TZ); 
  cmd.desiredThrustsN[1] = (1/(4.0)) * (forc_thrust - len_inv*TX + len_inv*TY + kVariable_inverse*TZ); 
  cmd.desiredThrustsN[2] = (1/(4.0)) * (forc_thrust + len_inv*TX - len_inv*TY + kVariable_inverse*TZ); 
  cmd.desiredThrustsN[3] = (1/(4.0)) * (forc_thrust - len_inv*TX - len_inv*TY - kVariable_inverse*TZ); 
  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  V3F momentCmd;
  const V3F diff_pqr = kpPQR * (pqrCmd - pqr);
  momentCmd = V3F(Ixx, Iyy, Izz) * diff_pqr;
  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  const float bXA = R(0,2), bYA = R(1,2);
  const float acc_thrust = -collThrustCmd / mass;
  const float bXC = accelCmd.x / (acc_thrust), bYC = accelCmd.y / (acc_thrust);
  const float diff_bXC = kpBank * (bXC - bXA);
  const float diff_bYC = kpBank * (bYC - bYA);
  const float inv_rot = 1.0F / R(2,2);
  pqrCmd.x = inv_rot * (R(1,0)*diff_bXC - R(0,0)*diff_bYC);
  pqrCmd.y = inv_rot * (R(1,1)*diff_bXC - R(0,1)*diff_bYC);
  pqrCmd.z = 0.0F;
  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;
  // Thrust 
  const float BZ = R(2,2);
  velZCmd = CONSTRAIN(velZCmd, -maxAscentRate, maxDescentRate);
  //err 
  const float err = posZCmd - posZ;
  const float diff_err = velZCmd - velZ;
  integratedAltitudeError += err * dt;
  const float U1_alter = kpPosZ * err + kpVelZ * diff_err + KiPosZ * integratedAltitudeError + accelZCmd;
  float acc_z_desired = (U1_alter - CONST_GRAVITY) / BZ;
  thrust = -acc_z_desired * mass;
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{ accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;
  V3F accelCmd = accelCmdFF;
  velCmd.x = CONSTRAIN(velCmd.x, -maxSpeedXY, maxSpeedXY);
  velCmd.y = CONSTRAIN(velCmd.y, -maxSpeedXY, maxSpeedXY);
  // Control Loops Params
  const V3F err = posCmd - pos;
  const V3F diff_err = velCmd - vel;
  accelCmd = kpPosXY*err + kpVelXY*diff_err + accelCmd;
  // Desired accel 
  accelCmd.x = CONSTRAIN(accelCmd.x, -maxAccelXY, maxAccelXY);
  accelCmd.y = CONSTRAIN(accelCmd.y, -maxAccelXY, maxAccelXY);
  accelCmd.z = 0.0F;
  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  float diff_yaw_CMd=0;
  const float err = normAng(yawCmd - yaw);
  diff_yaw_CMd = kpYaw * err;
  return diff_yaw_CMd;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
