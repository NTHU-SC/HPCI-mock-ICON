#! /usr/bin/perl
#
# header2installheader.pl.in --- apply config.h definitions on-disk
#
# Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
#
# Version: 1.0
# Author: Thomas Jahns <jahns@dkrz.de>
# Keywords:
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
use strict;
use warnings;

use File::Path qw(mkpath);
use File::Spec::Functions qw(catfile splitpath);
use File::Temp qw(tempfile);
use Getopt::Long ();
use cpp ();
use Storable qw(dclone);

sub checkPredicate(\%$$$);
sub substituteTokens($\@\%);
sub setup();
sub internConfigH();
sub applyInstDef();
sub hpp_transform($\@);


my ($debug, $verbose, $write_hpp_header,
    $srcBaseDir, $dstBaseDir,
    $parseCacheDB,
    $configHFn, %configHDefs, @configHInsert)
    = (exists($ENV{'DEBUG'})?int($ENV{'DEBUG'}):0, 0, 0,
       '', '');

our $headerLine;

my @symPrefixMacros = qw(SymPrefix SYMPREFIX symprefix);
my @symPrefixes = qw(PPM_ PPM_ ppm_);

setup();

internConfigH();

applyInstDef();

print(STDERR Data::Dumper->Dump([\%configHDefs], [qw(configHDefs)]))
    if ($debug > 0);
print(STDERR 'srcBaseDir = ', $srcBaseDir, "\n",
      'dstBaseDir = ', $dstBaseDir, "\n")
    if ($debug);

foreach my $headerFn (@ARGV)
{
  my $srcBaseDirStrLen = length($srcBaseDir);
  die('Unexpected prefix for ', $headerFn,
      ' should start with ', $srcBaseDir)
      if (substr($headerFn, 0, $srcBaseDirStrLen) ne $srcBaseDir);
  my ($instHeaderFn) = catfile($dstBaseDir,
                               substr($headerFn, $srcBaseDirStrLen+1));
  if ($write_hpp_header and substr($instHeaderFn, -4) eq '.hpp') {
    $headerFn = substr($headerFn, 0, -2);
  }
  my (undef, $instHeaderDir, undef) = splitpath($instHeaderFn);
  mkpath($instHeaderDir)
      unless (-d $instHeaderDir);
  die unless (-d $instHeaderDir);
  print(STDERR 'transforming ', $headerFn, ' to ', $instHeaderFn, "\n")
      if ($verbose > 0 or $debug > 0);
  my ($headerFh, $instHeaderFh);
  my $symprefixsubst = 0;
  open($headerFh, '<', $headerFn)
      or die('Cannot read header file: ', $headerFn, "\n");
  if (-e $instHeaderFn) {
    open($instHeaderFh, '+<', $instHeaderFn)
        or die('Cannot open output file for writing: ', $instHeaderFn, "\n");
  } else {
    open($instHeaderFh, '+>', $instHeaderFn)
        or die('Cannot open output file for writing: ', $instHeaderFn, "\n");
  }
  my @predicateLevels = ({'ignored' => 0, 'expand' => 0});
  my $ignored = 0;
  my (@headerContents) = (<$headerFh>);
  my (@headerLineNo) = (1..@headerContents);
  my $line = 0;
  while ($line < @headerContents)
  {
    local *headerLine = \$headerContents[$line];;
    my $lineno = $headerLineNo[$line];
    my $thisLineIgnored = $ignored;
    $headerLine = cpp::transliterateTrigraphs($headerLine);
    if ($headerLine
        =~ m{^\s*\#\s*(if|ifdef|ifndef|endif|else|elif)\b\s*(.*?)\s*$}xo)
    {
      # FIXME: handle multi-line conditionals
      my ($clause, $predicate) = ($1, $2);
      my ($predicateStartPos, $predicateEndPos) = ($-[2], $+[2]);
      print(STDERR 'Investigating ', $clause,
            (defined($predicate)?(' ', $predicate):()),
            ' at ', $headerFn, ', line ', $lineno, "\n") if $debug > 1;
      if ($clause eq 'ifdef' and $predicate eq 'HAVE_CONFIG_H')
      {
        my (@nextline);
        die('Unexpected #ifdef HAVE_CONFIG_H block structure at ', $headerFn,
            ', line ', $lineno, "\n")
            if (!defined($headerContents[$line+2])
                or !defined($headerContents[$line+1])
                or !$headerContents[$line+1]
                =~ m{^\#\s*include\s+(?:<|")config.h(?:"|>)}xo
                or !$headerContents[$line+2] =~ m{^\#\s*endif}xo);
        print(STDERR
              'Replacing #ifdef HAVE_CONFIG_H/#include "config.h"/#endif',
              ' structure at', "\n  ", $headerFn, ', line ', $lineno, "\n")
            if ($debug > 1);
        splice(@headerContents, $line, 3, @configHInsert);
        splice(@headerLineNo, $line, 3, (($lineno) x @configHInsert));
        redo;
      }
      print(STDERR 'line ', $lineno, "\n") if ($debug > 2);
      if ($clause eq 'endif')
      {
        die('#endif without opening #if at line ', $lineno,
            ' of ', $headerFn, "\n")
            if (@predicateLevels < 2);
        print(STDERR Data::Dumper->Dump([\@predicateLevels],
                                        ['predicatelevels']))
            if $debug > 1;
        print(STDERR '-' x 70, "\n",
              @headerContents[$line..($line+5>=@headerContents
                                      ?@headerContents-1:$line+5)],
              '-' x 70, "\n") if ($debug > 2);
        $thisLineIgnored = $predicateLevels[-1]{'expand'};
        pop(@predicateLevels);
        $ignored = $predicateLevels[-1]{'ignored'};
      }
      elsif ($clause eq 'else')
      {
        die('#else without opening #if at line ', $lineno,
            ' of ', $headerFn, "\n")
            if (@predicateLevels < 2);
        if (defined($predicateLevels[-1]{'ignored'}))
        {
          $ignored = $predicateLevels[-1]{'ignored'} =
              ($predicateLevels[-1]{'predValueAny'}
               or $predicateLevels[-2]{'ignored'});
        }
      }
      elsif ($clause eq 'elif')
      {
        die('#elif without opening #if at line ', $lineno,
            ' of ', $headerFn, "\n")
            if (@predicateLevels < 2);
        my ($predValue, $transformedPredicate)
            = checkPredicate(%configHDefs, $predicate, $headerFn, $lineno);
        if ((defined($predicateLevels[-1]{'predValueAny'})
             && $predicateLevels[-1]{'predValueAny'})
            || (defined($predValue) && !$predValue))
        {
          # This case is simple: the conditional can be precomputed
          # and is false or is in a block for which an expansion
          # already has been found. Consequently the conditional will
          # be removed and the block copied/ignored according to the
          # predicate evaluation
          $ignored = 1;
          $thisLineIgnored = 1;
          $predicateLevels[-1]{'ignored'} = $ignored;
        }
        elsif (!defined($predicateLevels[-1]{'predValueAny'})
               && defined($predValue) && $predValue)
        {
          # the previous conditionals were not computable and
          # this block is definitely true, that means we can transform
          # it into an else block and omit other alternatives
          $ignored = 0;
          $thisLineIgnored = 0;
          $predicateLevels[-1]{'predValueAny'} = 1;
          $headerLine
              =~ s{^\s*\#\s*(elif)\b\s*(.*?)\s*$}{\#else}xo;
        }
        elsif (defined($predicateLevels[-1]{'predValueAny'})
               && !$predicateLevels[-1]{'predValueAny'}
               && defined($predValue) && $predValue)
        {
          # the previous conditionals were all computable and false
          # and this block is definitely true, that means we can
          # transform it into an unconditional block
          $ignored = 0;
          $thisLineIgnored = 1;
          $predicateLevels[-1]{'predValueAny'} = 1;
          $predicateLevels[-1]{'expand'} = 1;
        }
        elsif (!defined($predicateLevels[-1]{'predValueAny'})
               && !defined($predValue))
        {
          # we don't know enough about this conditional, leave it in
          # place
          my $transformedHeaderLine
              = substr($headerLine, $predicateStartPos,
                       $predicateEndPos - $predicateStartPos,
                       $transformedPredicate);
          print(STDERR join("\n", 'Replacing', $headerLine, ' with',
                            $transformedHeaderLine), "\n")
              if $debug > 2;
          $headerLine = $transformedHeaderLine;
        }
        else
        {
          die('Unhandled conditional in ', $headerFn, ', line ', $lineno,
              ', this is a bug, please report to Thomas Jahns <jahns@dkrz.de>',
              "\n");
        }
      }
      else
      {
        if ($clause eq 'if') {
          # set expand to 1 if all symbols involved in $predicate are
          # defined/undefined in config.h
          my ($predValue, $transformedPredicate)
              = checkPredicate(%configHDefs, $predicate, $headerFn, $lineno);
          my $expand = defined($predValue)?1:0;
          $ignored = (($expand and !$predValue) or $ignored);
          push(@predicateLevels, { 'ignored' => $ignored,
                                   'predValue' => $predValue,
                                   'predValueAny' => $predValue,
                                   'line' => $lineno,
                                   'expand' => $expand,
                                 });
          $thisLineIgnored = $expand;
          print(STDERR 'In ', cpp::whoami(), ', line ', __LINE__, ': ',
                Data::Dumper->Dump([\@predicateLevels],
                                   ['predicatelevels']))
              if $debug > 2;
          if (!$thisLineIgnored) {
            # we don't know enough about this conditional, leave it in
            # place
            print(STDERR join("\n", 'Replacing predicate ', $predicate,
                              'in', $headerLine, ' with',
                              $transformedPredicate), "\n")
                if $debug > 2;
            substr($headerLine, $predicateStartPos,
                   $predicateEndPos - $predicateStartPos,
                   $transformedPredicate);
          }
        } elsif (($clause eq 'ifdef') or ($clause eq 'ifndef')) {
          die('Malformed predicate: ', $clause, ' ', $predicate)
              if (not cpp::validPreProcSym($predicate));
          my ($expand, $predValue) = (0, 1);
          if (exists($configHDefs{$predicate}))
          {
            $predValue = defined($configHDefs{$predicate});
            $predValue = !$predValue if ($clause eq 'ifndef');
            print(STDERR 'In ', cpp::whoami(), ', line ', __LINE__, ': ',
                  Data::Dumper->Dump([\$predValue], [qw(predicate)]))
                if $debug > 1;
            $ignored = (!$predValue or $ignored);
            $expand = 1;
            $thisLineIgnored = $expand;
          }
          push(@predicateLevels, { 'ignored' => $ignored,
                                   'predValue' => $predValue,
                                   'predValueAny' => $predValue,
                                   'line' => $lineno,
                                   'expand' => $expand});
          print(STDERR Data::Dumper->Dump([\@predicateLevels],
                                          ['predicatelevels']))
              if $debug > 2;
        }
      }
      print(STDERR ($ignored?'Ignoring':'Evaluating'), ' upcoming text.',
            "\n") if $debug > 1;
    } elsif (!$ignored) {
      if ($headerLine =~ m{^\s*#+(?:undef)\s+([A-Za-z_]\w*)\s*(..*)?$}o) {
        my ($preProcSym, $trail) = ($1, $2);
        die ($headerFn, ', line ', $lineno, ":\n",
             '#undef of config.h macro!', "\n")
            if (exists($configHDefs{$preProcSym}));
        warn('Trailing garbage on undef line ',
             $lineno, ' of ', $headerFn, ': ', $trail)
            if defined($trail);
      }
      # handle expansion of SymPrefix/SYMPREFIX/symprefix macros
      elsif ($headerLine =~ m{\#\s*include\s*[<"](?:core/)?symprefix.h[>"]}xo)
      {
        $thisLineIgnored = 1;
        $symprefixsubst = 1;
      }
      elsif ($headerLine =~ m{\#\s*include\s*([<"]).*?(?<=[^\\])\1}xo)
      {
        # keep other include directives unchanged
      }
      else
      {
        if ($symprefixsubst)
        {
          for(my $i =0; $i < @symPrefixMacros; ++$i)
          {
            $headerLine
                =~ s{$symPrefixMacros[$i]\((\w+)\)}{$symPrefixes[$i]$1}g;
          }
        }
        # tokenize and replace tokens from config.h definitions
        # FIXME: this simple approach will fail for function-like macros
        my @tokenMapping = eval { cpp::tokenize($headerLine); };
        if (!$@) {
          $headerLine = substituteTokens($headerLine, @tokenMapping,
                                         %configHDefs);
        }
      }
    }
    if ($ignored or $thisLineIgnored)
    {
      splice(@headerContents, $line, 1);
      splice(@headerLineNo, $line, 1);
    }
    else
    {
      ++$line;
    }
  }
  die('Unbalanced #if at line ', $predicateLevels[-1]{'line'},
      ' of ', $headerFn, "\n")
      if @predicateLevels > 2;
  if ($write_hpp_header) {
    hpp_transform($headerFn, @headerContents);
  }
  my (@oldHeaderContents) = (<$instHeaderFh>);
  unless (join('',@oldHeaderContents) eq join('', @headerContents)) {
    if ($debug and @oldHeaderContents) {
      no warnings 'once';
      local *STDTOUT;
      open(STDOUT, '>&', STDERR);
      if (open(my $diffFh, '|-', 'diff', '-u', $instHeaderFn, '-')) {
        print($diffFh @headerContents)
            or print (STDERR 'warning: problem writing to diff: ', $!, "\n");
        if (!close($diffFh) and $! != 0) {
          print (STDERR 'warning: problem running diff: ', $!, "\n");
        }
      } else {
        print (STDERR 'warning: could not launch diff: ', $!, "\n");
      }
      close(STDOUT);
    }
    truncate($instHeaderFh, 0)
        or die('cannot truncate ', $instHeaderFn, ": $!\n");
    seek($instHeaderFh, 0, 0)
        or die('failed to position output position: ', $!, "\n");
    print($instHeaderFh @headerContents);
  }
  close($headerFh) && close($instHeaderFh)
      or die;
}

sub internConfigH()
{
  my $configHFh;
  die('unknown path to config.h')
      unless (defined($configHFn) and -f $configHFn and -r $configHFn);
  open($configHFh, '<', $configHFn)
      or die('configuration header cannot be read', "\n");
  while (my $line = <$configHFh>)
  {
    if ($line =~ m{^/\*\s+\#\s*undef\s+(\w+)\s+\*/$}xo)
    {
      $configHDefs{$1} = undef;
    }
    elsif ($line =~ m{$cpp::definition_match}x)
    {
      my ($preProcSym, $args, $definition) = ($1, $2, $3);
      my ($lineno, $nextline) = ($.);
      while (substr($definition, -1) eq '\\'
             and $nextline = <$configHFh>)
      {
        $definition .= $nextline;
      }
      $configHDefs{$preProcSym} = cpp::parseDefinition($args, $definition,
                                                       $configHFn, $lineno);
    }
    elsif ($line =~ m{^\#\s*(?:if|ifdef|ifndef)\s+}xo)
    {
      # line is start of conditional block, copy it verbatim
      my @condBlock = ($line);
      my $currentCondLevel = 1;
      while ($line = <$configHFh>)
      {
        push(@condBlock, $line);
        if ($line =~ m{^\#\s*endif}xo)
        {
          last if (!--$currentCondLevel);
        }
        elsif ($line =~ m{^\#\s*(?:if|ifdef|ifndef)\s+}xo)
        {
          ++$currentCondLevel;
        }
      }
      die('incomplete conditional block', "\n")
          if (!defined($line));
      push(@configHInsert, @condBlock);
    }
    else
    {
      print(STDERR 'Ignoring line: ', $line)
          if ($debug > 1);
    }
  }
}

# test whether result of preprocessor conditional can be determined in advance
sub checkPredicate(\%$$$)
{
  my ($defines, $predicate, $headerFn, $lineno) = @_;
  my @tokenMapping = cpp::tokenize($predicate);
  print(STDERR Data::Dumper->Dump([\@tokenMapping],
                                  ['tokenMapping']))
      if $debug > 2;
  my $result = 4711;
  foreach my $token (@{$tokenMapping[0]})
  {
    print(STDERR 'Inspecting token ', $token, ' at ', __LINE__+1, "\n")
        if $debug > 4;
    next if ($token eq 'defined');
    print(STDERR 'Inspecting token ', $token, ' at ', __LINE__+1, "\n")
        if $debug > 4;
    next if (exists($defines->{$token}));
    print(STDERR 'Inspecting token ', $token, ' at ', __LINE__+1, "\n")
        if $debug > 4;
    next if ($token =~ m{$cpp::operator_match});
    print(STDERR 'Inspecting token ', $token, ' at ', __LINE__+1, "\n")
        if $debug > 4;
    next if ($token =~ m{$cpp::integral_constant_match});
    print(STDERR 'This predicate cannot be deduced because of token "',
          $token, '"', "\n")
        if $debug > 3;
    $result = undef;
    last;
  }
  if (defined($result))
  {
    print(STDERR 'This predicate could be deduced', "\n")
        if $debug > 2;
    $result = eval {
      cpp::evaluatePreProcExpression($predicate, $defines);
    };
    if ($@)
    {
      die($headerFn, ':', $lineno, ":\n",
          'Error: invalid conditional expression: ',
          $predicate, "\n");
    }
    return ($result);
  } else {
    return (undef, substituteTokens($predicate, @tokenMapping, %$defines));
  }
}

sub substituteTokens($\@\%)
{
  my ($headerLine, $tokenMapping, $configHDefs) = @_;
  # FIXME: when tokenizing doesn't work that usually means a multi-line
  # construct needs to be accounted for
  print(STDERR 'In ', cpp::whoami(), ', line ', __LINE__, ' ',
        Data::Dumper->Dump([$tokenMapping], ['tokenMapping']))
      if $debug > 3;
  # list of replacements, pre-filled with original tokens
  $tokenMapping->[2] = [@{$tokenMapping->[0]}];
  # list of adjusted replacement positions, pre-filled with original
  # token positions
  $tokenMapping->[3] = dclone($tokenMapping->[1]);
  # flag if replaced, initially 0 and then 1 for every replacement
  $tokenMapping->[4] = [ (0) x @{$tokenMapping->[1]} ];
  for(my $tok = 0; $tok < @{$tokenMapping->[0]}; ++$tok) {
    my $token = $tokenMapping->[0][$tok];
    if (defined($token) && exists($configHDefs->{$token})
            && defined($configHDefs->{$token}))            {
      print STDERR 'Substituting ',
          join(' ', @{$configHDefs->{$token}{'definition'}[0]}), ' for ',
          $token, "\n" if $debug > 1;
      $tokenMapping->[4][$tok] = 1;
      my $replacement = join(' ', @{$configHDefs->{$token}{'definition'}[0]});
      $tokenMapping->[2][$tok] = $replacement;
      my $lenChange = length($replacement)
          - length($tokenMapping->[0][$tok]);
      substr($headerLine, $tokenMapping->[3][$tok][0],
             $tokenMapping->[3][$tok][1] + 1 - $tokenMapping->[3][$tok][0])
          =~ s/$token/@{$configHDefs->{$token}{'definition'}[0]}/;
      for (my $tokRest = $tok + 1;
           $tokRest < @{$tokenMapping->[0]};
           ++$tokRest)
      {
        $tokenMapping->[3][$tokRest][0] += $lenChange;
        $tokenMapping->[3][$tokRest][1] += $lenChange;
        print(STDERR Data::Dumper->Dump([$tokenMapping],
                                        ['tokenMapping']))
            if $debug > 3;
      }
    }
  }
  return $headerLine;
}

# and now for something completely different: find function
# declarations that contain C99 VLA pointers and fix them for C++
# which insists that C89/C90 is the only C there is.  Let's not think
# too hard how C++ could just support something C's had for 25 years
# now (and is a widespread extension already in the early 90ies).
sub hpp_transform($\@)
{
  my ($headerFn, $headerContents, $numLines, $insert_start, $insert_end, $fundecl_seen, $restrict_seen)
      = ($_[0], $_[1], scalar(@{$_[1]}), -1, -1, 0, 0);
  my (undef, undef, $headerBaseFn) = File::Spec->splitpath($headerFn);
  print STDERR $headerBaseFn, "\n" if ($debug > 2);
  my $fname_stanza_re = qr([@\\]file\s+${headerBaseFn});
  for (my $i = 0;
       $i < $numLines;
       ++$i)
  {
    my ($line) = ($headerContents->[$i]);
    if ($line =~ m{^typedef\b}o) {
      # handle typedef block
      print STDERR 'typedef', "\n" if ($debug > 1);
    } elsif ($line =~ m{^\#\s*include\s*\<xt/.*\.h\>}xo) {
      print STDERR 'internal include', "\n" if ($debug > 1);
      print(STDERR $line, $+[0], "\n") if ($debug > 2);
      substr($headerContents->[$i], $+[0]-1, 0) = 'pp';
    } elsif (index($line, '/*') != -1) {
      print STDERR 'comment', $i, "\n" if ($debug > 1);
      # begin of comment, scan to end of comment
      while ($line !~ m{\*/}o) {
        $line = $headerContents->[++$i];
        if ($line =~ m{$fname_stanza_re}) {
          print STDERR 'doxygen filename stanza', "\n" if ($debug > 1);
          substr($headerContents->[$i], $+[0], 0) = 'pp';
        }
      }
      print STDERR 'endcomment', $i, "\n" if ($debug > 1);
      print STDERR Data::Dumper->Dump([\$line], [qw(line)])
          if ($debug > 2);
    } elsif ($line =~ m{^([a-zA-Z_][0-9a-zA-Z_]*\s+)
                        *(?:[a-zA-Z_][0-9a-zA-Z_]*)(?=\()}gxo)
    {
      $fundecl_seen = 1;
      print STDERR 'fundecl ', $line,
          ' at line ', $i, "\n" if ($debug > 3);
      # start of function declaration found, identify
      # range up to its end
      my $preceding_context = 0;
      while ($i > $preceding_context
             and index($headerContents->[$i - $preceding_context - 1], ';') == -1
             and index($headerContents->[$i - $preceding_context - 1], '*/') == -1)
      {
        ++$preceding_context;
      }
      # identify function declarations with C99 array arguments and mangle
      my ($fn_start, $fn_end) = ($i-$preceding_context);
      if ($insert_start == -1) {
        $insert_start = $fn_start;
      }
      while ($i < @$headerContents
             and index($headerContents->[$i], ';') == -1) { ++$i; }
      $insert_end = $fn_end = $i;
      $line = join('', @{$headerContents}[$fn_start..$fn_end]);
      $line =~ m{^(?:[a-zA-Z_][0-9a-zA-Z_]*\s+)*([a-zA-Z_][0-9a-zA-Z_]*)(?=\()}gmo;
      print STDERR 'fundecl ', $1,
          ($debug > 2 ? (' at lines ', $fn_start+1, ' to ', $fn_end) : ()),
          "\n" if ($debug > 1);
      my @params;
      while ($line =~ m{\G([(,]\s*+([^]()[,]+
                            (?:\[[a-z_A-Z0-9]*\]|\(\)|(?1)\))*)
                          \s*(?=[),]))}gxo) {
        print STDERR "match: '", $&, "'\n" if ($debug > 2);
        push(@params, [ $-[2], $+[2], $2]);
      }
      print STDERR Data::Dumper->Dump([\$line, \@params], [qw(line params)])
          if ($debug > 2);
      for (my $i = 0;
           $i < @params;
           ++$i) {
        $params[$i][1] = $params[$i][1] - $params[$i][0];
      }
      my ($vla_re, $vla_dim_name_re, $needs_hpp_replacement) = (undef, 0);
      foreach my $param (@params) {
        if ((' '.$param->[2].' ') =~ m{\brestrict\b}o) {
          $restrict_seen = 1;
          print(STDERR 'restrict seen in param "', $param->[2], "\"\n")
              if ($debug > 1);
        }
        if ($param->[2] ne 'void'
            and $param->[2] =~ m{^(?:[A-Za-z_][A-Za-z0-9_]*
                                 \s+)+([A-Za-z_][A-Za-z0-9_]*)$}xo) {
          $param->[3] = $1;
          if (defined($vla_dim_name_re)) {
            $vla_dim_name_re = qr{$vla_dim_name_re|$1};
          } else {
            $vla_dim_name_re = qr{$1};
          }
          $vla_re = qr{\s+\**([A-Za-z_][A-Za-z0-9_]*)
                       (?:\[\])*\[($vla_dim_name_re)(?:\[.*\])*\]}x;
        } elsif (defined($vla_re)
                 and index($param->[2], '[') != -1
                 and $param->[2] =~ m{$vla_re}x) {
          $param->[4] = join('', substr($param->[2], 0, $-[1]),
                             '*', $1);
          $needs_hpp_replacement = 1;
        }
      }
      if ($needs_hpp_replacement) {
        my $ofs = 0;
        foreach my $param (@params) {
          if (defined($param->[4])) {
            my $ofs_adj = length($param->[4]) - $param->[1];
            substr($line, $param->[0]+$ofs, $param->[1]) = $param->[4];
            $ofs += $ofs_adj;
          }
        }
      }
      print STDERR Data::Dumper->Dump([\$line, \@params], [qw(line params)])
          if ($debug > 2);
      splice(@$headerContents, $fn_start, $fn_end - $fn_start + 1,
             split(/^/m, $line));
    }
  }
  if ($fundecl_seen) {
    my @extra_hdr_contents;
    @extra_hdr_contents = ('#ifndef restrict', '#define restrict', '#endif')
        if ($restrict_seen);
    push(@extra_hdr_contents, 'extern "C" {');
    splice(@$headerContents, $insert_end+1, 0, '};' . "\n");
    splice(@$headerContents, $insert_start, 0,
           map({ $_ . "\n" }
               @extra_hdr_contents));
  }
}


{

  # choose default appropriate for LP64
  my ($c_sizeof_int, $c_sizeof_long, $c_sizeof_long_long,
      $c_char_is_unsigned) = (4, 8, 8, 0);

  sub setup()
  {
    parseOptions();
    cpp::init({ 'c_sizeof_int' => $c_sizeof_int,
                'c_sizeof_long' => $c_sizeof_long,
                'c_sizeof_long_long' => $c_sizeof_long_long,
                'c_char_is_unsigned'
                => $c_char_is_unsigned,
              });
  }

  sub parseOptions
  {
    local $_;
    my ($help, $usage)=(0, 0);
    my $optionParser = new Getopt::Long::Parser;
    # This program should accept the full set of preprocessor flags,
    # but might not implement all of them. Therefore abbreviation and
    # ignoring of case are not allowed because either might promote
    # unhandled options to options of this program. Also unhandled
    # options are passed through.
    configureOptionParser($optionParser, 'no_auto_abbrev', 'no_ignore_case');
    # duplicate -- option terminator for brain-dead versions of
    # Getopt::Long (which unfortunately includes the one distributed
    # with Perl 5.6.1)
    if (grep /^--$/, @ARGV)
    {
      my @seen = (0, 0);
      @ARGV = ((grep { ($seen[0] or /^--$/)?
                           ($seen[0]?0:!($seen[0]=1)):1 } @ARGV),
               '--', '--',
               (grep { ($seen[1] or /^--$/)?$seen[1]++:0 } @ARGV));
    }
    my $result
        = $optionParser->getoptions('debug+' => \$debug,
                                    'help|?!' => \$help,
                                    'usage!' => \$usage,
                                    'verbose+' => \$verbose,
                                    'parse-cache=s' => \$parseCacheDB,
                                    'config-header=s' => \$configHFn,
                                    'srcdir=s' => \$srcBaseDir,
                                    'dstdir=s' => \$dstBaseDir,
                                    'symbol-prefix-c=s' => \$symPrefixes[0],
                                    'symbol-prefix-F=s' => \$symPrefixes[1],
                                    'symbol-prefix-f=s' => \$symPrefixes[2],
                                    'inst-def|D=s' => \&set_inst_def,
                                    'inst-undef|U=s' => \&set_inst_def,
                                    'c-cast-mode=s' => \&set_cast_mode,
                                    'c-sizeof-int=i' => \$c_sizeof_int,
                                    'c-sizeof-long=i' => \$c_sizeof_long,
                                    'c-sizeof-long-long=i'
                                      => \$c_sizeof_long_long,
                                    'c-char-is-signed'
                                    => \&set_c_char_sign,
                                    'c-char-is-unsigned'
                                    => \&set_c_char_sign,
                                    'write-hpp-header!'
                                    => \$write_hpp_header
                                   );
    if($help or $usage)
    {
      pod2usage( {
                  '-msg' => "",
                  '-exitval' => 0,
                  '-verbose' => 1
                 });
    }
    $cpp::debug = $debug;
  }

   sub set_cast_mode($)
   {
     my $mode = $_[1];
     if ($mode =~ m{^lp64$}ixo)
     {
       ($c_sizeof_int, $c_sizeof_long, $c_sizeof_long_long) = (4, 8, 8);
     }
     elsif ($mode =~ m{^ilp32$}ixo)
     {
       ($c_sizeof_int, $c_sizeof_long, $c_sizeof_long_long) = (4, 4, 8);
     }
     else
     {
       print STDERR 'Unexpected arithmetic conversion mode: ', $mode, "\n";
       exit(1);
     }
   }

   my (%extraInstDefs, %extraInstUndefs);

   sub set_inst_def($)
   {
     my ($opt_name, $opt_val) = @_;
     if ($opt_name eq 'inst-def')
     {
       my ($macro, $definition) = ($opt_val =~ m{([^=]*)(?:=(.*))?}o);
       if (defined($macro))
       {
         $definition = '' unless (defined($definition));
         $extraInstDefs{$macro} = $definition;
       }
       else
       {
         die('Unexpected installation define ', $opt_val, "\n");
       }
     }
     elsif ($opt_name eq 'inst-undef')
     {
       die('Unexpected installation define ', $opt_val, "\n")
           if (index($opt_val, '=') != -1);
       $extraInstUndefs{$opt_val} = undef;
     }
   }

   sub applyInstDef()
   {
     while (my ($macro, $definition) = each(%extraInstDefs))
     {
       $configHDefs{$macro} = $definition;
     }
     foreach my $macro (keys(%extraInstUndefs))
     {
       delete($configHDefs{$macro});
     }
   }

   sub set_c_char_sign($$)
   {
      my $mode = $_[0];
      if ($mode eq 'c-char-is-unsigned')
      {
         $c_char_is_unsigned = 1;
      }
      elsif ($mode eq 'c-char-is-signed')
      {
         $c_char_is_unsigned = 0;
      }
      else
      {
         die('Unexpected error');
      }
   }

}

# evil hack to work with older versions of Getopt::Long
sub configureOptionParser
{
   my ($optionParser, @options) = (@_);
   eval {
      $optionParser->configure(@options);
   };
   if ($@)
   {
      my $save = Getopt::Long::Configure ($optionParser->{settings}, @options);
      $optionParser->{settings} = Getopt::Long::Configure($save);
   }
}


# Local Variables:
# mode: cperl
# cperl-indent-level: 2
# license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:
