warning('off','MATLAB:mex:GccVersion_link')

% store the path from which this function is called
current_path = pwd;

% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);


fprintf('Compiling LOBPCG\n');

error_msg = 'The C compiler could not succesfully compile ';

flags = '-largeArrayDims';
cflags = 'CFLAGS="\$CFLAGS -std=c99 -O3 -DLAPACK"';
flags = [cflags ' ' flags];

s = sprintf ('mex %s -O -lm -lmwlapack %s.c', flags, 'lobpcg_mex') ;
eval (s);
% if mex('-outdir', current_path, fullfile(current_path,'lobpcg_mex.c'),...
%          '-lm', '-lmwlapack' , '-largeArrayDims',...
%         'CFLAGS="\$CFLAGS -std=c99 -g"')
%     error([error_msg, mex_path]);
% end