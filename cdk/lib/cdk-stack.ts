import * as ecr from '@aws-cdk/aws-ecr';
import * as cdk from '@aws-cdk/core';
import { core } from '@myhelix/cdk-library';
import { Namer } from "multi-convention-namer";
import { Batch } from './batch';
import { Local } from './local';

export class PGXStack extends core.Stack {
  constructor(scope: cdk.Construct, id: string, props: core.StackProps) {
    super(scope, id, props);

    const local = new Local();

    // s3 bucket
    const bucket = new Namer([
        this.account,
        'hipaa-exome-workflow',
        PGXStack.nameModifier(this.namedEnv),
    ]).kebab;

    // ecr
    const ecrArn = `arn:aws:ecr:${this.region}:${local.ecrAccountId}:repository/helix-py-app-r2v`;
    const ecrRepository = ecr.Repository.fromRepositoryArn(this, 'ecrRepo', ecrArn);
    const pgxImageTag = process.env.CIRCLE_SHA1 || 'latest';
    const pgxVersion = process.env.CIRCLE_TAG || 'unversioned';

    // batch
    const batch = new Batch(this, 'Batch', {
      vpc: this.vpc,
      bucket: bucket,
      ecrRepository: ecrRepository,
      projectName: local.projectName,
      imageTag: pgxImageTag,
      r2vVersion: pgxVersion,
      environment: this.namedEnv,
      route53record: local.route53Record
    });

  }

  static nameModifier(namedEnv: core.NamedEnv): string {
    // deploying on the master account requires an additional qualifier to the default names
    return namedEnv.account == core.Environment.MasterProduction.account
      ? namedEnv.name
      : "";
  }
}
