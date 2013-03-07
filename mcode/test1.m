function test1
   figure('KeyPressFcn',@printfig);
   function printfig(src,evnt)
      if evnt.Character == 'e'
         print ('-deps',['-f' num2str(src)])
      elseif length(evnt.Modifier) == 1 & strcmp(evnt.Modifier{:},'control') & evnt.Key == 't'
         print ('-dtiff','-r200',['-f' num2str(src)])
   end

