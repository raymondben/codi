function [codi_score, check_score]=codi(rtimes, sao, verbose);
% function [codi_score, check_score]=codi(rtimes, sao, verbose);
%
% The CODI algorithm.
%
% Inputs:
%  rtimes - times of occurrence of RR intervals, in seconds
%  sao - the pulse oximetry saturation signal, sampled at 1 Hz
%  verbose=0 - 0 [default, no output],1 [print final summary only], or 2 [also show progress]
%
% Returns:
%  codi_score - the CODI score
%  check_score - the difference between the autonomic arousal index and the resaturation index. Large negative
%                values indicate many more arousals than resaturations; large positive values the opposite.

% Ben Raymond

% check inputs
if nargin<3
  verbose=0;
end
if nargin<2
  error('two inputs expected');
end

% calculate inter-beat interval
rri=diff(rtimes); % RR interval
rtimes=rtimes(2:end); % (ending) time of each interval

rec_length=(rtimes(end)-rtimes(1))/3600;

if min(size(rtimes))>1,
  error('R-times data must be a vector.')
elseif min(size(rri))>1,
  error('RRI data must be a vector.')
end

% some algorithm parameters
dip_delay=40; % max number of seconds delay between cardiac arousal and sao dip
rriaro_tachy_thresh=-0.06; % threshold for tachycardia (absolute change in RR interval, in seconds)
rriaro_brady_thresh=0.085; % threshold for bradycardia (absolute change in RR interval, in seconds)
rriaro_winlen=10; % length of window to use for finding tachy/bradycardias
rriaro_max_delay=20; % maximum delay allowed between tachy and brady to count as arousal (in seconds)
rriaro_min_sep=10; % minimum separation between successive arousals for them to be considered separate (in seconds)

if verbose>1, fprintf('Evenly resampling RRI signal ... ');end
rri=spline(rtimes,rri,0:0.1:rtimes(end)); % resample at 10Hz rate, starting from t=0
rri=rri(1:10:end); % downsample to 1Hz
if verbose>1, fprintf('done.\n');end

if verbose>1, fprintf('Finding cardiac arousals ... ');end
rriaro=fastrriaro(rri,rriaro_tachy_thresh,rriaro_brady_thresh,rriaro_winlen,rriaro_max_delay);
if verbose>1, fprintf('done.\n');end

% delete arousals not separated by at least rriaro_min_sep
rriaro(:,find(diff(rriaro(1,:))<rriaro_min_sep)+1)=[];
arousal_times=rriaro(1,:);
arousal_times=arousal_times-1; % because we resampled to 1Hz, starting at t=0
    
% match cardiac arousals with sao disturbances
if verbose>1, fprintf('Correlating cardiac and oximetry events ... ');end  
%find rises in sao signal
[risecount,saorise_times]=saocodi(sao);
% find cardiac arousals which are followed by dip
codi_times=[];
codi_score=0;
for l=1:length(arousal_times),
  temp=saorise_times-arousal_times(l);
  temp=temp(find(temp>=0)); 
  temp=temp(find(temp<=dip_delay));
  if ~isempty(temp),
    codi_score=codi_score+1;
    if ~isempty(codi_times),
      codi_times=[codi_times arousal_times(l)];
    else
      codi_times=arousal_times(l);
    end
  end
end
% convert codi_score to hourly index
codi_score=codi_score/rec_length;
if verbose>1, fprintf('done.\n');end

% check difference between arousals and resaturations
check_score=(length(arousal_times)-length(saorise_times))/rec_length;
        
%calculate dip counts for Sa02 signal
if verbose>1, fprintf('Calculating Sa02 dip counts ... ');end
dips=dipcount_strad(sao); % use stradling's algorithm
if verbose>1, fprintf('done.\n\n\n');end

if verbose, 
  % display some results
  fprintf('CODI analysis (%s)\n----\n',datestr(now));
  fprintf('Recording length: %.1f hours\n',rec_length);
  fprintf('CODI score: %.1f\n',codi_score);
  fprintf('Difference between autonomic arousal index and resaturation index: %.1f\n',check_score);
  fprintf('----\nOximetry-only analysis\n');
  fprintf('4%% dip index: %.1f /hour\n',dips(4));
  fprintf('3%% dip index: %.1f /hour\n',dips(3));
  fprintf('2%% dip index: %.1f /hour\n',dips(2));
  fprintf('----\n\n');
end
    

  
  
  
  
function f=fastrriaro(rri,tachy_thresh,brady_thresh,winlen,max_delay);
% f=function fastrriaro(rri,tachy_thresh,brady_thresh,winlen,max_delay);
%
% An attempt to find RRI arousals in a timely fashion!
%
% Inputs:
%   rri - RRI signal, evenly resampled at 1 Hz
%   tachy_thresh - threshold for tachycardia, in seconds (abs change in RRI)
%   brady_thresh - threshold for bradycardia, in seconds
%   winlen - length of window to use for finding tachy/bradycardias
%   max_delay - max delay allowed between tachy and brady to count as arousal
%
%  Returns:
%    2xn matrix, 1st row is the sample number of arousal, 2nd row is the magnitude of change in the RR interval during the arousal period

verbose=1;
if nargin<1
  error('Not enough input arguments');
elseif nargin<5,
  if verbose, 
    fprintf('fastrriaro: using default parameters\n');
  end
  max_delay=20; % max allowable delay between end of tachy and start of brady to class as arousal
  tachy_thresh=-0.06;
  brady_thresh=0.085;
  winlen=10;
end
neg_sl_thresh=-1*abs(tachy_thresh)/winlen;
pos_sl_thresh=abs(brady_thresh)/winlen;

graphic_disp=0;
debug_on=0;

% make sure rri is a (row) vector
rrisize=size(rri);
if min(rrisize)>1,
  error('rri input must be a vector')
end
% assume rri sampled at 1 Hz
rrilen=max(rrisize);
if rrisize(2)==1,
  rri=rri';
end

% prefilter rri
%[h3,g3,rh3,rg3]=daub(4);
% OR use 
h3=[-0.12940952255126   0.22414386804201   0.83651630373781   0.48296291314453];
rh3=[0.48296291314453   0.83651630373781   0.22414386804201  -0.12940952255126];
rrf=aprox(rri,h3,rh3,3); %use 3-scales approximation to rri

% find slope of filtered signal
sl=slidingslope(rrf,winlen);
% actually want sl to reflect the slope of the trailing window
% adjust offset of sl to approximate this
% sl will be a row variable
sl=[sl(1)*ones(1,floor(winlen/2)) sl(1:length(sl)-floor(winlen/2))];

neg_sl=find(sl<=neg_sl_thresh);

% retain starts and ends of runs only
pd=diff(neg_sl);
endel=[]; startel=[];
pdo=find(pd==1);
pdosave=pdo;
if pdo(1)==1,
  pdo(1)=[];
  startel=[1];
end
if ~isempty(startel),
  startel=[startel pdo(find(pd(pdo-1)~=1))];
else
  startel=[pdo(find(pd(pdo-1)~=1))];
end
pdo=pdosave;
clear pdosave
if pdo(length(pdo))==length(pd), 
  endel=pdo(length(pdo))+1;
  pdo(length(pdo))=[];
end
endel=[pdo(find(pd(pdo+1)~=1))+1 endel];

% endel and startel ought to be the same length but it is not guaranteed
neg_sl=[neg_sl(startel);neg_sl(endel)];

clear pd pdo endel temp_neg_sl


pos_sl=find(sl>=pos_sl_thresh);
% find start and end of runs
nd=diff(pos_sl);
startel=[]; endel=[];
ndo=find(nd==1); 
ndosave=ndo;
if ndo(1)==1,
  ndo(1)=[];
  startel=[1];
end
startel=[startel ndo(find(nd(ndo-1)~=1))];
ndo=ndosave;
clear ndosave
if ndo(length(ndo))==length(nd), 
  endel=ndo(length(ndo))+1;
  ndo(length(ndo))=[];
end
if ~isempty(endel),
  endel=[ndo(find(nd(ndo+1)~=1))+1 endel];
else
  endel=[ndo(find(nd(ndo+1)~=1))+1];
end

pos_sl=[pos_sl(startel);pos_sl(endel)];

clear nd ndo startel temp_pos_sl

% correlate these two sets of slopes
f=[-1 -1]'; % temp data only there
for k=1:size(neg_sl,2),
  temp=pos_sl(1,:)-neg_sl(2,k);
  tempidx=(find(temp>0));
  temp=temp(tempidx);
  if ~isempty(temp),
    posind=tempidx(find(temp<=max_delay));
    if ~isempty(posind)
      posind=posind(1);
      if isempty(find(abs(f(1,:)-(neg_sl(1,k)))<1e-06)),
        % find size of tachycardia
        if verbose==2,
          fprintf('negative slope from %g to %g, pos from %g to %g\n',neg_sl(1,k),neg_sl(2,k),pos_sl(1,posind),pos_sl(2,posind))
        end
        maxtrrf=max(rrf(neg_sl(1,k):pos_sl(1,posind))); % max during tachy
        mintrrf=min(rrf(neg_sl(1,k):pos_sl(1,posind))); % min during tachy 
        f=[f [neg_sl(1,k);maxtrrf-mintrrf]];
      end
    end
  end
end

f(:,1)=[];
%f(1,:)=f(1,:)-winlen; %change so that start of arousal (approx) is returned
return  




function [saocount,saoepoch]=saocodi(sao,samp_rate);
% function [count,times]=saorise(sao,find_resat_only,samp_rate,return_raw_data);
% 
% Finds the count of rises in the Sa02 signal. Used by the CODI script.
%  
% Inputs:
%   sao - Sa02 signal (range 0-100, in units of %)
%   samp_rate - sampling rate (default=1Hz)
%
% Returns:
%   count - number of rises per hour
%   times - vector containing list of times (in samples) at which a 
%           significant rise in Sa02 was found
  
% Ben Raymond, November 2000

% returns times in samples not epochs

  verbose=0;

  % check inputs
  if nargin<2, 
    samp_rate=1;
  end
  
  if ~isempty(find(sao>100)) | ~isempty(find(sao<0))
    error('SaO2 signal has elements beyond 100 or below 0.')
  end
  
  % define parameters
  saosl_thresh=0.15;
  % define window width to use for finding rise rate
  wind=10;
  if verbose, fprintf('\n saorise_s.m: threshold %g over %g s\n',saosl_thresh*wind,wind),end
  
  if size(sao,1)==1,
    sao=sao';
  end
  
  % find slope of sao
  saosl=slidingslope(sao,wind);
  % want to threshold this data
  exthresh=find(saosl>=saosl_thresh);
  %  need to remove consecutive runs here
  saoepoch=run_ends(exthresh);
  if ~isempty(saoepoch),
    % keep first of runs only
    saoepoch=exthresh(saoepoch(1,:));
    saocount=length(saoepoch);    
    % now convert to hourly rate
    saocount=saocount/length(sao)*samp_rate*3600;
  else
    saocount=0;
  end  
  return
  

  
  

function f=dipcount_strad(sao,samp_rate);
% function dc=dipcount_strad(sao,samp_rate);
% 
% Finds the dip count of Sa02 signal based on John Stradling's algorithm
%  
% Inputs:
%   sao - Sa02 signal
%   samp_rate - sampling rate (default=1Hz)
%
% Returns:
%   f - number of 1,2,3 and 4% dips per hour (row vector)
    
    
% Ben Raymond, June 1999
  
  if nargin<2, samp_rate=1;end

  if samp_rate~=1,
    error('Sample rates other than 1 Hz not implemented!')
  end
  
  debug=0;
  verbose=1;
  
  % convert from fraction to percent if needed
  if mean(sao)<1, %is fraction not percent
    if verbose, 
      fprintf('dipcount.m: Converting Sa02 to percent values ... ')
    end
    sao=100*sao;
    if verbose, fprintf('done.\n'),end
  end

  % first resample at 12s period, keep min in each period
  % NOTE ASSUMES SAMPLE RATE OF 1 HZ !!
  newsao=[];
  for k=1:12:length(sao)-11,
    newsao=[newsao min(sao(k:k+11))];
  end
  sao=newsao;
    
  dip_thresh=[1 2 3 4];
  rise_thresh=[1 1 2 3];
  f=[];
  
  for run_num=1:size(dip_thresh,2),
    dip_count=0;
    current_max=-1;
    k=1;
    while k<=length(sao),
      current_max=max(current_max,sao(k));
      % find dips of >dip_thresh%
      if sao(k)<(current_max-dip_thresh(run_num)),
  	dip_count=dip_count+1;
  	current_min=sao(k);
  	% find resat of >rise_thresh%
  	while k<=length(sao),
	  if sao(k)>current_min+rise_thresh(run_num),
	    current_max=sao(k);
	    break
	  end
	  k=k+1;
  	end
      end
      k=k+1;
    end
    if ~isempty(f),
      f=[f dip_count/length(sao)/12*3600];
    else
      f=[dip_count/length(sao)/12*3600];
    end
  end
  return
    
  
function f=slidingslope(y,window);
% function f=slidingslope(y,window);
%
% Finds the slope of vector y over a sliding window of length window
% (uses best least-squares fit)
% Assumes that y lies in the middle of the window (best to use odd-length window)
% Much faster than using nlfilter or similar

% Ben Raymond, November 2000

% check inputs
if min(size(y))>1, error('Input data must be a vector');end

if window>length(y), error('Window length too large');end

if window==1, f=y;return;end

% make sure y is a col vector
if size(y,1)==1, 
  y=y';
  ywasrow=1;
else
  ywasrow=0;
end

% make matrix Y, cols are sliding observations of y
firstcol=[1:window]';
firstrow=[0:length(y)-window];

idx=firstcol(:,ones(1,length(firstrow))) + firstrow(ones(window,1),:);
Y=y(idx);

% demean columns of Y
Y=Y-repmat(sum(Y)/window,window,1);

% construct matching X matrix and demean
X=[1:window]'; 
X=(X-mean(X)); 
X=repmat(X,1,size(Y,2));


f=(sum(X .*Y) ./ sum(X.^2));

% add initial and ending slopes to make lengths match
f=[f(1)*ones(1,floor(window/2)) f f(length(f))*ones(1,window-floor(window/2)-1)];

if ~ywasrow,f=f';end
return  



function f=run_ends(x);
% function f=run_ends(x);
%
% Returns the element numbers which are the first and last in runs of consecutive integers
%
% Inputs:
%   x - the data to check
%
% Returns:
%   f - row 1 is the starting elements of runs; row 2 is the corresponding ending elements

% Ben Raymond, June 2000

x=reshape(x,1,length(x));

pd=diff(x);
endel=[]; startel=[];
pdo=find(pd==1);
if ~isempty(pdo)
  pdosave=pdo;
  if pdo(1)==1,
    pdo(1)=[];
    startel=[1];
  end
  startel=[startel pdo(find(pd(pdo-1)~=1))];
  pdo=pdosave;
  clear pdosave
  if pdo(length(pdo))==length(pd), 
    endel=pdo(length(pdo))+1;
    pdo(length(pdo))=[];
  end
  endel=[pdo(find(pd(pdo+1)~=1))+1 endel];
  
  f=[startel;endel];
else
  f=[];
end
return

  
  
function y=aprox(x,h,rh,sc);

%  APROX    Obtain a projection on the approximation space in the context
%           of the multiresolution analysis. 
%	
%	    Y = APROX (X,H,RH,SC) calculates the approximation of signal 
%           in X at the scale SC, using the analysis lowpass filter H and
%           the synthesis one RH.
%
%	    If APj is the approximation of X at scale 2^j, and DEi 
%           are the DETAIL outputs, for i=1...j, then
%
%		X = AP  + DE + ... + DE
%		      j	    j          1
%
%	    See also: DETAIL, MULTIRES, MRES2D, WT.	


%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
%                                                      
%                                                      
% Uvi_Wave is free software; you can redistribute it and/or modify it      
% under the terms of the GNU General Public License as published by the    
% Free Software Foundation; either version 2, or (at your option) any      
% later version.                                                           
%                                                                          
% Uvi_Wave is distributed in the hope that it will be useful, but WITHOUT  
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    
% for more details.                                                        
%                                                                          
% You should have received a copy of the GNU General Public License        
% along with Uvi_Wave; see the file COPYING.  If not, write to the Free    
% Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.             
%                                                                          
%       Author: Sergio J. Garcia Galan
%               Nuria Gonzalez Prelcic
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------

x=x(:)';
h=h(:)';
rh=rh(:)';

lx=length(x);
lf=length(h)+length(rh);
dd=floor((lf)/2)-1;
d=dd;
if sc>1
  for i=2:sc
    d=d+dd*2^(i-1);
  end
end

tt=x;

for i=1:sc,
  tt=conv(h,tt);
  tt=tt(1:2:length(tt));
end
for i=1:sc,
  tt=[tt;zeros(1,length(tt))];
  tt=conv(rh,tt(:)');
end

y=tt(1+d:lx+d);
return
