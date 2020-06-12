% function mtx_to_tbl_pgfplots(filename,comment,legend_str,data,[precision])
%
% Save the data in a matrix to a text file which can be read by pgfplots.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename: the data is saved in a file with the name filename.table
%    size: string
% comment: A single line comment about the data
%    size: string
% legend_str: a string with a space separated list of legends pertaining to
% the data columns
%    size: string
% data: a data matrix
%    size: NxM
% precision: an optional scalar or vector of the column precision used when
% saving the data. The default value is a precision of 5 for each column.
%    size: 1x1 or 1xM or Mx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% void
%
% v. 1.1.0
% Created:        2013-12-18  by Jesper Kjær Nielsen (jkn@es.aau.dk)
% Last modified:  2013-12-19  by Jesper Kjær Nielsen (jkn@es.aau.dk)
function mtx_to_tbl_pgfplots(filename,comment,legend_str,data,precision)
  % data matrix size
  [N,M] = size(data);
  % print a warning if N is very large
  if N > 1000
    warning('There are more than 1000 rows in the data matrix. LaTeX is not going to be happy. Please consider downsampling your data matrix.');
  end
  % set the precision if not set
  if nargin < 5
    % if precision is not set, set it to 5
    precision = 5*ones(1,M);
  elseif length(precision) == 1
    % use the same precision in all columns
    precision = precision*ones(1,M);
  end
  % open a file for writing
  fid = fopen([filename,'.table'], 'wt');
  % add the comment to the file
  fprintf(fid,['%%', comment,'\n']);
  % add the legends to the file
  fprintf(fid,[legend_str,'\n']);
  % write the data to the file
  for n=1:N
    for m=1:M
      fprintf(fid,'%.*f ',precision(m),data(n,m));
    end
    fprintf(fid,'\n');
  end
  fclose(fid);
end
