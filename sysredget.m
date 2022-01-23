function o = sysredget(options,name,default,flag)
%SYSREDGET Get SYSRED OPTIONS parameters.
%   VAL = SYSREDGET(OPTIONS,'NAME') extracts the value of the named parameter
%   from SYSRED options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%   
%   VAL = SYSREDGET(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified (is [])
%   in OPTIONS.  For example
%   
%     val = SYSREDget(opts,'Tolred',1e-4);
%   
%   returns val = 1e-4 if the Tolred property is not specified in opts.
%   
%   See also SYSREDSET.

%   RELEASE 2.0 of SLICOT Model and Controller Reduction Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Revisions:
%   V. Sima, Research Institute for Informatics, Bucharest, Feb. 2009.

if nargin < 2
   error('Not enough input arguments.');
end
if nargin < 3
   default = [];
end
if nargin < 4
   flag = [];
end

% undocumented usage for fast access with no error checking
if isequal('fast',flag)
   o = sysredgetfast(options,name,default);
   return
end

if ~isempty(options) && ~isa(options,'struct')
   error('First argument must be an options structure created with SYSREDSET.');
end

if isempty(options)
   o = default;
   return;
end

optionsstruct = sysredset;

Names = fieldnames(optionsstruct);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error(['Unrecognized property name ''%s''.  ' ...
         'See SYSREDSET for possibilities.'], name);
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' Names{j(1),:}];
    for k = j(2:length(j))'
      msg = [msg ', ' Names{k,:}];
    end
    error('%s).', msg);
  end
end

if any(strcmp(Names,Names{j,:}))
  o = getfield(options, Names{j,:});
  if isempty(o)
    o = default;
  end
else
  o = default;
end

%------------------------------------------------------------------
function value = sysredgetfast(options,name,defaultopt)
%SYSREDGETFAST Get SYSRED OPTIONS parameter with no error checking so fast.
%   VAL = SYSREDGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is 
%   probably a subset of the options in OPTIONS.
%

S.type = '.'; S.subs = name;

if ~isempty(options)
    value = subsref(options,S);
else
    value = [];
end

if isempty(value)
    value = subsref(defaultopt,S);
end


%-------------------------------------------------
function f = getfield(s,field)
%GETFIELD Get structure field contents.
%   F = GETFIELD(S,'field') returns the contents of the specified
%   field.  This is equivalent to the syntax F = S.field.
%   S must be a 1-by-1 structure and field must be deblanked already.  
% 

sref.type = '.'; sref.subs = field;
f = subsref(s,sref);
