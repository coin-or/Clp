@REM This file will run command line tests, comparable to test\makefile.am

@echo off
SET "BINDIR=%~1"
SET "SAMPLEDIR=%~2"
SET "NETLIBDIR=%~3"
SET "MIPLIBDIR=%~4"

echo INFO: Running %0 %*
echo INFO: Using bindir '%BINDIR%' and sampledir '%SAMPLEDIR%' and netlibdir '%NETLIBDIR%'
echo INFO: Miplibdir '%MIPLIBDIR%' is ignored for these tests.

if "%BINDIR%"=="" echo ERROR: No bindir given. && goto :usage
if not exist "%BINDIR%" echo ERROR: Folder bindir %BINDIR% does not exist. && goto :usage

if "%SAMPLEDIR%"=="" echo ERROR: No sampledir given. && goto :usage
if not exist %SAMPLEDIR% echo ERROR: Folder sampledir %SAMPLEDIR% does not exist. && goto :usage
if not %errorlevel%==0 echo ERROR: %SAMPLEDIR% cannot contain spaces. && goto :usage

if "%NETLIBDIR%"=="" echo ERROR: No netlibdir given. && goto :usage
if not exist %NETLIBDIR% echo ERROR: Folder netlibdir %NETLIBDIR% does not exist. && goto :usage
if not %errorlevel%==0 echo ERROR: %NETLIBDIR% cannot contain spaces. && goto :usage

goto :test

:usage
echo INFO: Usage %0 ^<bindir^> ^<sampledir^> ^<netlibdir^>
echo INFO: where ^<bindir^> contains the executables, sampledir the sample files and netlibdir the netlib files.
echo INFO: This script runs automated test.
echo INFO: The ^<sampledir^> and ^<netlibdir^> must not contain spaces!
echo INFO: For example: %0 "D:\Some Directory\" ..\..\samples\ "..\..\netlib"
goto :error

:error
echo ERROR: An error occurred while running the tests!
echo INFO: Finished Tests but failed
exit /b 1

:test
echo INFO: Starting Tests

"%BINDIR%\clp.exe" -dirSample %SAMPLEDIR% -unitTest -dirNetlib %netlibdir% -netlib
if not %errorlevel%==0 echo ERROR: Error running clp.exe tests. && goto :error

"%BINDIR%\osiUnitTest.exe" -mpsDir=%SAMPLEDIR% -netlibDir=%NETLIBDIR% -testOsiSolverInterface 
if not %errorlevel%==0 echo ERROR: Error running osiUnitTest.exe tests. && goto :error

echo INFO: Finished Tests successfully (%ERRORLEVEL%)