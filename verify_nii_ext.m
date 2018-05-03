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

%  Verify NIFTI header extension to make sure that each extension section
%  must be an integer multiple of 16 byte long that includes the first 8
%  bytes of esize and ecode. If the length of extension section is not the
%  above mentioned case, edata should be padded with all 0.
%
%  Usage: [ext, esize_total] = verify_nii_ext(ext)
%
%  ext - Structure of NIFTI header extension, which includes num_ext,
%       and all the extended header sections in the header extension.
%       Each extended header section will have its esize, ecode, and
%       edata, where edata can be plain text, xml, or any raw data
%       that was saved in the extended header section.
%
%  esize_total - Sum of all esize variable in all header sections.
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function [ext, esize_total] = verify_nii_ext(ext)

   if ~isfield(ext, 'section')
      error('Incorrect NIFTI header extension structure.');
   elseif ~isfield(ext, 'num_ext')
      ext.num_ext = length(ext.section);
   elseif ~isfield(ext, 'extension')
      ext.extension = [1 0 0 0];
   end

   esize_total = 0;

   for i=1:ext.num_ext
      if ~isfield(ext.section(i), 'ecode') | ~isfield(ext.section(i), 'edata')
         error('Incorrect NIFTI header extension structure.');
      end

      ext.section(i).esize = ceil((length(ext.section(i).edata)+8)/16)*16;
      ext.section(i).edata = ...
	[ext.section(i).edata ...
	 zeros(1,ext.section(i).esize-length(ext.section(i).edata)-8)];
      esize_total = esize_total + ext.section(i).esize;
   end

   return                                       % verify_nii_ext

