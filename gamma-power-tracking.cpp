/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

#include "gamma-power-tracking.h"
#include <iostream>
#include <main_window.h>
#include <math.h>

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new GammaPowerTracking();
}

static DefaultGUIModel::variable_t vars[] = {
  // { "GUI label", "Tooltip description",
  //   DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  { "gamma-LFP", "gamma-LFP (V)", DefaultGUIModel::OUTPUT, },
  { "LFP", "LFP (A)", DefaultGUIModel::INPUT, },
  { "Wavelet1 freq. (Hz)", "Wavelet1 freq. (Hz)",
    DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
  { "Wavelet2 freq. (Hz)", "Wavelet2 freq. (Hz)",
    DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
  {"A State", "Tooltip description", DefaultGUIModel::STATE, },
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

GammaPowerTracking::GammaPowerTracking(void)
  : DefaultGUIModel("GammaPowerTracking with Custom GUI", ::vars, ::num_vars)
{
  initParameters();
  setWhatsThis("<p><b>GammaPowerTracking:</b><br> [...] </p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  update(INIT); 
  refresh();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

GammaPowerTracking::~GammaPowerTracking(void)
{
}

void
GammaPowerTracking::execute(void)
{
  /* EVERYTHING SHOULD BE TRANSLATED T0 SI UNITS HERE */    
  systime = count * period; // time in seconds for display

  double real = 0;
  double imag = 0;
  double wvl_norm = 0;
  mean_LFP_new = 0;
  for(int i=0;i<Length_wavelet1-1;i++)
    {
      LFP_history_vector[i+1] = LFP_history_vector[i]; //move all element to the left except first one
      increase_Real_and_Imaginary_of_Morlet_TF(&real, &imag,
					       LFP_history_vector[i+1]-mean_LFP_last,
					       100., 1e-4,  w0,
					       i, Length_wavelet1);
      mean_LFP_new += LFP_history_vector[i]/Length_wavelet1;
    }
  LFP_history_vector[0] = input(0); // update the history vector with the new recorded value
  mean_LFP_new += LFP_history_vector[0]/Length_wavelet1;
  mean_LFP_last = mean_LFP_new;
  output(0) = mean_LFP_last;
  count++;
}

void
GammaPowerTracking::increase_Real_and_Imaginary_of_Morlet_TF(double* real, double* imag, double X,
								    double freq, double dt, double w0,
								    int it, int N_wavelet_vector)
{
  double factor = 0;
  factor = (w0/2./sqrt(2.*M_PI)/freq)*(1.+exp(-pow(w0,2)/2)) * exp(-0.5 * (pow(((2.*M_PI*freq*(it-N_wavelet_vector)*dt)/w0),2)));
  *real += X * cos(2. * M_PI * freq * (it-N_wavelet_vector) * dt) * factor;
  *imag += - X * sin(2. * M_PI * freq * (it-N_wavelet_vector) * dt) * factor;
}



void
GammaPowerTracking::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-9; // time in s
      setParameter("Wavelet1 freq. (Hz)", wavelet1_freq);
      setParameter("Wavelet2 freq. (Hz)", wavelet2_freq);
      // setParameter("Wavelet3 freq. (Hz)", wavelet3_freq);
      initWavelets();
      // setState("A State", Length_wavelet1);
      break;

    case MODIFY:
      wavelet1_freq = getParameter("Wavelet1 freq. (Hz)").toDouble();
      wavelet2_freq = getParameter("Wavelet2 freq. (Hz)").toDouble();
      // wavelet3_freq = getParameter("Wavelet3 freq. (Hz)").toDouble();
      initWavelets();
      break;

    case UNPAUSE:
      break;

    case PAUSE:
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-9; // time in s
      break;

    default:
      break;
  }
}

void
GammaPowerTracking::initParameters(void)
{
  period = RT::System::getInstance()->getPeriod() * 1e-9; // time in s
  count = 0;
  systime = 0;
  wavelet1_freq = 100.;
  wavelet2_freq = 70.;
}

void
GammaPowerTracking::initWavelets(void)
{
  Length_wavelet1 = static_cast<int> (ceil(2 * pow(2, .5) * (w0/(M_PI*wavelet1_freq)) / period));
  std::cout << Length_wavelet1 << '\n';
  std::cout << wavelet1_freq << '\n';
}
