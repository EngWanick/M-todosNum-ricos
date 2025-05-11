@echo off

REM === Caminhos ===
SET PYTHON_EXE=C:\Program Files\Python312\python.exe
SET SIMULADOR=D:\Engenharia Mecânica\Mestrado\2º Sem\MetNum\M-todosNum-ricos\Programa_1\sedimentacao.py
SET GRAFICOS=D:\Engenharia Mecânica\Mestrado\2º Sem\MetNum\M-todosNum-ricos\Programa_1\gerar_graficos.py

REM === Executa a simulação ===
echo Iniciando simulação numérica...
call "%PYTHON_EXE%" "%SIMULADOR%"
echo.
echo ====> Simulação concluída. CSV gerado com sucesso.
echo.

REM === Pergunta ao usuário se deseja gerar os gráficos ===
set /p gerar="Deseja gerar os gráficos agora? (s/n): "
if /I "%gerar%"=="s" (
    echo Executando script de geração de gráficos...
    call "%PYTHON_EXE%" "%GRAFICOS%"
    echo Gráficos gerados.
) else (
    echo Operação finalizada sem geração de gráficos.
)

pause
