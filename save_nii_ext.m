%Copyright (c) 2014, Jimmy Shen
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

%  Save NIFTI header extension.
%
%  Usage: save_nii_ext(ext, fid)
%
%  ext - struct with NIFTI header extension fields.
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function save_nii_ext(ext, fid)

   if ~exist('ext','var') | ~exist('fid','var')
      error('Usage: save_nii_ext(ext, fid)');
   end

   if ~isfield(ext,'extension') | ~isfield(ext,'section') | ~isfield(ext,'num_ext')
      error('Wrong header extension');
   end

   write_ext(ext, fid);

   return;                                      % save_nii_ext


%---------------------------------------------------------------------
function write_ext(ext, fid)

   fwrite(fid, ext.extension, 'uchar');

   for i=1:ext.num_ext
      fwrite(fid, ext.section(i).esize, 'int32');
      fwrite(fid, ext.section(i).ecode, 'int32');
      fwrite(fid, ext.section(i).edata, 'uchar');
   end

   return;                                      % write_ext

