@echo off
echo Iniciando simulação de sedimentação...

set PYTHON=python
set SCRIPT=sedimentacao_fluido_particula_lucaswanick_242104541.py

if not exist %SCRIPT% (
    echo Erro: arquivo %SCRIPT% não encontrado!
    exit /b 1
)

%PYTHON% %SCRIPT%

if exist stokes_curves_final.csv (
    echo Arquivo CSV gerado com sucesso!
) else (
    echo Erro: CSV não foi gerado.
    exit /b 2
)

set /p RESPOSTA=Deseja visualizar os graficos agora? [s/n] 
if "%RESPOSTA%"=="s" (
    %PYTHON% visualizar_graficos.py
)

echo Execução finalizada.
pause
