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



%  Load NIFTI header extension after its header is loaded using load_nii_hdr.
%
%  Usage: ext = load_nii_ext(filename)
%
%  filename - NIFTI file name.
%
%  Returned values:
%
%  ext - Structure of NIFTI header extension, which includes num_ext,
%       and all the extended header sections in the header extension.
%       Each extended header section will have its esize, ecode, and
%       edata, where edata can be plain text, xml, or any raw data
%       that was saved in the extended header section.
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function ext = load_nii_ext(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: ext = load_nii_ext(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
   vox_offset = 0;
    
   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         if new_ext
            fseek(fid,108,'bof');
            vox_offset = fread(fid,1,'float32');
         end

         ext = read_extension(fid, vox_offset);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               error(msg);
            end

            if new_ext
               fseek(fid,108,'bof');
               vox_offset = fread(fid,1,'float32');
            end

            ext = read_extension(fid, vox_offset);
            fclose(fid);
         end
      end
   end

   return                                       % load_nii_ext


%---------------------------------------------------------------------
function ext = read_extension(fid, vox_offset)

   ext = [];

   if vox_offset
      end_of_ext = vox_offset;
   else
      fseek(fid, 0, 'eof');
      end_of_ext = ftell(fid);
   end

   if end_of_ext > 352
      fseek(fid, 348, 'bof');
      ext.extension = fread(fid,4)';
   end

   if isempty(ext) | ext.extension(1) == 0
      ext = [];
      return;
   end

   i = 1;

   while(ftell(fid) < end_of_ext)
      ext.section(i).esize = fread(fid,1,'int32');
      ext.section(i).ecode = fread(fid,1,'int32');
      ext.section(i).edata = char(fread(fid,ext.section(i).esize-8)');
      i = i + 1;
   end

   ext.num_ext = length(ext.section);

   return                                               % read_extension

