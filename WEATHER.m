%% WEATHER.m
% ***********************************************************************
% This file is part of the uibkCARNOT Blockset.
% 
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%   Dietmar Siegele     dietmar.siegele@uibk.ac.at
%   Eleonora Leonardi   eleonora.leonardi@uibk.ac.at
% Additional Copyright for this file see list auf authors.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its 
%    contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
% THE POSSIBILITY OF SUCH DAMAGE.
% **********************************************************************
% 
%% carnotUIBK version 1.3
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS,EL     2017-03-12   initial revision v1.0

%%
classdef WEATHER
    % WEATHER
    
    properties
        name = '';                  % name of the weather
        path = '';                  % path to weather file
        time_value = 0;
        latitude = 0;               % 1 latitude [°]   
        longitude = 0;              % 2 longitude [°]
        latitude_timezone = 0;      % 3 latitude difference of position to latitude time zone  [°]
        radiation_beam_normal = 0;  % 4 radiation beam normal [W/m²]
        radiation_diffuse_horizontal = 0;  % 5 radiation diffuse horizontal [W/m²]
        t_ambient = 0;              % 6 Ambient temperature [°C]
        t_sky = 0;                  % 7 Sky Temperature [°C]
        rh = 0;                     % 8 relative humidity [%]
        precip = 0;                 % 9 precipitation [m/s]
        cloud = 0;                  % 10 cloud index from 0 to 1
        p = 0;                      % 11 pressure [Pa]
        vw = 0;                     % 12 windspeed [m/s]
        wdir = 0;                   % 13 wind direction 0 = north, 90 = east
        incidence = 0;              % 14 incidence angle sun [°]
        tetap = 0;                  % 15 incidence angle longitudinal [°]
        tetas = 0;                  % 16 incidence angle trasversal [°]
        Idirect_surface = 0;        % 17 solar radiation direct surface [W/m²]
        Idiffuse_surface = 0;       % 18 solar radiation diffuse surface [W/m²]
        zenith = 0;                 % 19 zenith [°]   
        azimuth = 0;                % 20 azimuth [°] 
        
    end
    
    methods
        function obj = WEATHER(name, time_value, zenith, azimuth, latitude, longitude, latitude_timezone, radiation_beam_normal, radiation_diffuse_horizontal, t_ambient, t_sky, rh, precip, cloud, p, vw, wdir, incidence, tetap, tetas, Idirect_surface, Idiffuse_surface, path)
            obj.name = name;
            obj.time_value = time_value;
            obj.zenith = zenith;
            obj.azimuth = azimuth;
            obj.latitude = latitude;
            obj.longitude = longitude;
            obj.latitude_timezone = latitude_timezone;
            obj.radiation_beam_normal = radiation_beam_normal;
            obj.radiation_diffuse_horizontal = radiation_diffuse_horizontal;
            obj.t_ambient = t_ambient;
            obj.t_sky = t_sky;
            obj.rh = rh;
            obj.precip = precip;
            obj.cloud = cloud;
            obj.p = p;
            obj.vw = vw;
            obj.wdir = wdir;
            obj.incidence = incidence;
            obj.tetap = tetap;
            obj.tetas = tetas;
            obj.Idirect_surface = Idirect_surface;
            obj.Idiffuse_surface = Idiffuse_surface;
            obj.path = path;
            
        end
        
        function plot(obj)
            % to plot the external temperature, the beam normal radiation
            % and the diffusive radiation of the climate
            figure
            name = ['WEATHER' obj.name ': ambient temperature'];
            we = timeseries(obj.t_ambient(:,2), obj.t_ambient(:,1), 'Name', name);
            we.TimeInfo.StartDate = '01-01-2014';
            we.TimeInfo.Format = 'dd-mmm';
            we.DataInfo.Units = '°C';
            plot(we)
            xlabel('Time')
            ylabel('Temperature (°C)')
            title(name)
            grid on
            xlim([datetime('31-Dec-2016 23:59:00') datetime('01-Jan-2018 00:00:00')])
            
            figure
            name = ['WEATHER' obj.name ': radiation beam normal'];
            we = timeseries(obj.radiation_beam_normal(:,2), obj.radiation_beam_normal(:,1), 'Name', name);
            we.TimeInfo.StartDate = '01-01-2014';
            we.TimeInfo.Format = 'dd-mmm';
            we.DataInfo.Units = 'W/m²';
            plot(we)
            xlabel('Time')
            ylabel('Radiation (W/m²)')
            title(name)
            grid on
            xlim([datetime('31-Dec-2016 23:59:00') datetime('01-Jan-2018 00:00:00')])

            figure
            name = ['WEATHER' obj.name ': diffuse radiation horizontal'];
            we = timeseries(obj.radiation_diffuse_horizontal(:,2), obj.radiation_diffuse_horizontal(:,1), 'Name', name);
            we.TimeInfo.StartDate = '01-01-2014';
            we.TimeInfo.Format = 'dd-mmm';
            we.DataInfo.Units = 'W/m²';
            plot(we)
            xlabel('Time')
            ylabel('Radiation (W/m²)')
            title(name)
            grid on
            xlim([datetime('31-Dec-2016 23:59:00') datetime('01-Jan-2018 00:00:00')])
        end

        function weat_mon = calculate_weather(obj)
            T_m = 365*24*3600 + ([0 31 (28+31) (31+28+31) (30+31+28+31) (31+30+31+28+31) (30+31+30+31+28+31) (31+30+31+30+31+28+31) (31+31+30+31+30+31+28+31) (30+31+31+30+31+30+31+28+31) (31+30+31+31+30+31+30+31+28+31) (30+31+30+31+31+30+31+30+31+28+31) (31+30+31+30+31+31+30+31+30+31+28+31)]*24*3600);
            T_m = T_m/(obj.time_value(2,1)-obj.time_value(1,1));
            for jj = 1:12
                weat_mon(jj) = mean(obj.t_ambient((T_m(jj)+1):(T_m(jj+1)),2));
            end
        end
    end
end