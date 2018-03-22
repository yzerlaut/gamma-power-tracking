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
  { "pLFP", "pLFP (V)", DefaultGUIModel::OUTPUT, },
  { "LFP", "LFP (A)", DefaultGUIModel::INPUT, },
  { "Wavelet freq. (Hz)", "Wavelet freq. (Hz)",
    DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
  { "Smoothing (ms)", "Smoothing (ms)",
    DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
  // { "Wavelet2 freq. (Hz)", "Wavelet2 freq. (Hz)",
  //   DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
  {"A State", "Tooltip description", DefaultGUIModel::STATE, },
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

GammaPowerTracking::GammaPowerTracking(void)
  : DefaultGUIModel("GammaPowerTracking", ::vars, ::num_vars)
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
  mean_LFP_new = 0;

  for(int i=Length_wavelet1-1;i>0;i--)
    {
      LFP_history_vector[i] = LFP_history_vector[i-1]; //move all element to the left except first one
      increase_Real_and_Imaginary_of_Morlet_TF(&real, &imag,
					       LFP_history_vector[i]-mean_LFP_last,
					       wavelet1_freq, period,  w0, wavelet1_norm,
					       i, Length_wavelet1);
      mean_LFP_new += LFP_history_vector[i]/Length_wavelet1;
    }

  LFP_history_vector[0] = input(0); // update the history vector with the new recorded value
  mean_LFP_new += LFP_history_vector[0]/Length_wavelet1;
  mean_LFP_last = mean_LFP_new;
  cum_pLFP += sqrt(pow(real,2)+pow(imag,2))/ismoothing;

  if (count%ismoothing==0) 
    { // every smoothing update
      // we start the output calculus
      output(0) = 0 ;
      for(int i=Length_pLFP_vector-1;i>0;i--) 
	{
	  pLFP_history_vector[i] = pLFP_history_vector[i-1];
	  output(0) += pLFP_history_vector[i]*pLFP_history_norm_vector[i] ;
	}

      pLFP_history_vector[0] = cum_pLFP; 
      output(0) += pLFP_history_vector[0]*pLFP_history_norm_vector[0];
      cum_pLFP = 0; // we reset pLFP

    }
  count++;
}

void
GammaPowerTracking::increase_Real_and_Imaginary_of_Morlet_TF(double* real, double* imag, double X,
							     double freq, double dt, double w0, double norm_factor,
							     int it, int N_wavelet_vector)
{
  double factor = exp(-0.5 * (pow(((2.*M_PI*freq*(it-N_wavelet_vector/2)*dt)/w0),2)));
  *real += X * cos(2. * M_PI * freq * (it-N_wavelet_vector/2) * dt) * factor / norm_factor *dt;
  *imag += - X * sin(2. * M_PI * freq * (it-N_wavelet_vector/2) * dt) * factor / norm_factor *dt;
}

void
GammaPowerTracking::normalization_factor(double* norm_factor,
					 double freq, double dt, double w0)
{
  *norm_factor = (w0/2./sqrt(2.*M_PI)/freq)*(1.+exp(-pow(w0,2)/2));
}



void
GammaPowerTracking::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      setParameter("Wavelet freq. (Hz)", wavelet1_freq);
      setParameter("Smoothing (ms)", Tsmoothing);
      // std::cout << period << '\n'; /!\ Be careful with the acquisisiton frequency, in the POSIX mode, 1kHz /!\
      // setParameter("Wavelet2 freq. (Hz)", wavelet2_freq);
      // setParameter("Wavelet3 freq. (Hz)", wavelet3_freq);
      // setState("A State", Length_wavelet1);
      period = RT::System::getInstance()->getPeriod() * 1e-9; // time in s
      initWavelets();
      break;

    case MODIFY:
      period = RT::System::getInstance()->getPeriod() * 1e-9; // time in s
      wavelet1_freq = getParameter("Wavelet freq. (Hz)").toDouble();
      Tsmoothing =  getParameter("Smoothing (ms)").toDouble();
      // wavelet2_freq = getParameter("Wavelet2 freq. (Hz)").toDouble();
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
  count = 0;
  systime = 0;
  Tsmoothing = 20.;
  wavelet1_freq = 100.;
  cum_pLFP = 0; // we reset pLFP
}

void
GammaPowerTracking::initWavelets(void)
{
  Length_wavelet1 = static_cast<int> (ceil(2 * sqrt(2) * (w0/M_PI/wavelet1_freq) / period));
  normalization_factor(&wavelet1_norm, wavelet1_freq, period, w0);
  ismoothing = static_cast<int> (ceil(Tsmoothing * 1e-3 / period / 2)); // Tsmoothing in ms !
  cum_pLFP = 0; // we reset pLFP
  // finding the exponential normalization
  double full_norm = 0;
  for (int i=0;i<Length_pLFP_vector;i++) {
    pLFP_history_norm_vector[i] = exp(-i/(2*ismoothing));
    full_norm += pLFP_history_norm_vector[i];
      }
  for (int i=0;i<Length_pLFP_vector;i++) {
    pLFP_history_norm_vector[i] /= full_norm;
  }
  // std::cout << Length_wavelet1 << '\n';
  // std::cout << Length_wavelet1 ;
}
