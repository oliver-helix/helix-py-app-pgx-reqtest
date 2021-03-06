#!/usr/bin/env node
import * as cdk from '@aws-cdk/core';
import { core } from "@myhelix/cdk-library";
import { PGXStack } from '../lib/pgx-stack';
import { Local } from '../lib/local';
import { Namer } from 'multi-convention-namer';

const app = new cdk.App();
const local = new Local();

// derive environment
const namedEnv = core.Environment.findByName(process.env.ENVIRONMENT);

const stackName = new Namer([local.projectName, namedEnv.account == core.Environment.MasterProduction.account ? namedEnv.name : '', 'stack']).camel;

new PGXStack(
    app,
    stackName,
    {
        accountingTag: core.AccountingCategory.ENGINEERING,
        serviceTag: `${namedEnv.name}-${local.projectName}-cdk`,
        namedEnv: namedEnv
    }
);
