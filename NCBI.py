import ut_functions

def getApiKey():
    params = ut_functions.readIni("pipeline.ini")
    return (params['api_key'])

def fixStatement(statement, syst=""):
    statement = statement.replace("esearch ", "esearch -api %s " % getApiKey())
    statement = statement.replace("efetch ", "efetch -api %s " % getApiKey())
    return (statement)